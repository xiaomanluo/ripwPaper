library("tidyverse")
library("survival")
library("lfe")
library("survminer")
source("ripw.R")

## Preprocess data
data_full <- read.csv("../data/opentable_us_misc.csv") %>%
        filter(state.name != "District of Columbia") %>%
        select(state.name, date,
               reserv.diff,death, confirmed, Population, 
               Neighbor_soe.status, Neighbor_tot,
               soe.status, BEDS, dem_to_rep_ratio) %>%
        rename(treat = soe.status,
               state = state.name,
               vote = dem_to_rep_ratio,
               beds = BEDS / Population)%>%
        mutate(date = as.Date(date),
               confirmed = log(confirmed + 1),
               beds = log(beds)) %>%
        select(-contains("Neighbor")) %>%
        na.omit()

start_date <- as.Date("2020-02-29")
end_date <- as.Date("2020-03-13")
T <- as.integer(end_date - start_date) + 1
data <- data_full %>%
    mutate(time = date - start_date) %>%
    filter(date >= start_date & date <= end_date)

## cross-fitting
set.seed(2021)
nfolds <- 10
nstates <- length(unique(data$state))
states <- unique(data$state)
foldid <- gen_cf_inds(nstates, nfolds)

## Generate data for coxph
survdata <- list()
for (i in 1:length(states)){
    tmpdata <- data %>%
        filter(state == states[i]) %>%
        select(treat, time, confirmed, beds, vote)
    inds <- which(tmpdata$treat > 0)
    if (length(inds) > 0){
        ind <- min(inds)
        tmpdata <- tmpdata[1:ind, ] %>%
            mutate(time1 = time - 1,
                   time2 = time) %>%
            mutate(status = 0,
                   id = i)
        tmpdata$status[length(tmpdata$status)] <- 1
    } else {
        tmpdata <- tmpdata %>%
            mutate(time1 = time - 1,
                   time2 = time) %>%
            mutate(status = 0,
                   id = i)
    }
    survdata[[i]] <- tmpdata
}
survdata <- do.call(rbind, survdata) %>%
    mutate(time1 = as.numeric(time1),
           time2 = as.numeric(time2),
           time = as.numeric(time)) %>%
    select(confirmed, beds, vote, time, time1, time2, status, id) 

## test for proportional hazard assumption
fit <- coxph(Surv(time1, time2, status) ~ confirmed + vote + beds, data = survdata)
test <- cox.zph(fit)
ggcoxzph(test)

# prediction by state
gps_est <- rep(NA, nstates)
adopt_time <- rep(NA, nstates)
adopt <- rep(NA, nstates)

for (k in 1:nfolds){
    ids <- foldid[[k]]
    train <- filter(survdata, !(id %in% ids))
    test <- filter(survdata, id %in% ids)
    status_df <- test %>% group_by(id) %>%
        summarize(status = max(status))
    adopt[ids] <- status_df$status
    fit <- coxph(Surv(time1, time2, status) ~ confirmed + vote + beds, data = train)
    obj <- summary(survfit(fit, newdata = test, id = id))
    for (index in ids){
        status <- status_df$status[status_df$id == index]
        if (is.null(obj$strata)){
            survtime <- obj$time
            survprob <- obj$surv
        } else {
            survtime <- obj$time[obj$strata == index]
            survprob <- obj$surv[obj$strata == index]
        }
        if (length(survprob) == 0){
            ## If adopt_time[index] is earlier than
            ## adoption times of all training units,
            ## use P(T <= start_date) instead
            earliest_date <- min(summary(survfit(fit))$time)
            newtest <- test %>% filter(id == index)
            added_row <- newtest[nrow(newtest), ]
            added_row$time2 <- earliest_date
            added_row$time1 <- earliest_date - 1
            added_row$status <- 1
            newtest <- rbind(newtest, added_row)
            newobj <- summary(survfit(fit, newdata = newtest, id = id))
            survprob <- newobj$surv
            survtime <- T
        }
        if (status == 1){
            if (length(survprob) == 1){
                gps_est[index] <- 1 - survprob
            } else {
                gps_est[index] <- -diff(tail(survprob, 2))
            }
        } else {
            ## For always control units, we regularize
            ## the gps as P(T >= end_date - 1)
            gps_est[index] <- tail(survprob, 2)[1]
        }
        if (status){
            adopt_time[index] <- max(survtime)
        }
    }
}

## Reshaped distribution
adopt_time[is.na(adopt_time)] <- T + 1
adopt_time <- T + 1 - adopt_time

Pi <- rep(0, T + 1)
Pi[2:T] <- 1 / 2 / T
Pi[c(1, T+1)] <- (T + 1) / 4 / T

EPi <- head(cumsum(rev(Pi)), -1)
Theta <- Pi[adopt_time + 1] / gps_est

## Generate the outcomes and treatment assignments
Y <- data %>% select(state, time, reserv.diff) %>%
    spread(time, reserv.diff) %>%
    select(-state) %>%
    as.matrix
tr <- data %>% select(state, time, treat) %>%
    spread(time, treat) %>%
    select(-state) %>%
    as.matrix

## Unweighted two-way estimator     
obj_unw <- felm(reserv.diff ~ treat + confirmed + vote | factor(state) + factor(time) | 0 | 0,
                data = data,
                exactDOF = TRUE,
                cmethod = "reghdfe")
unw_tauhat <- obj_unw$coefficients[1]
unw_se <- summary(obj_unw)$coefficients[1,2]
unw_pval <- 2 * pnorm(abs(unw_tauhat / unw_se), lower.tail = FALSE)

## RIPW estimator
X <- data %>%
    select(treat, confirmed, vote) %>%
    as.matrix
X <- lapply(1:T, function(j){
    X[(j-1)*nstates + 1:nstates, ]
})
X <- Reduce(cbind, X)
muhat <- twfe_mhat(Y, tr, X,
                   timevarying = TRUE,
                   joint = TRUE,
                   foldid = foldid)
obj_ripw <- ripw(Y, tr, muhat = muhat$mhat0, Theta = Theta)

## Summary
res <- data.frame(
    unw_tauhat = unw_tauhat, unw_se = unw_se, unw_pval = unw_pval,
    ripw_tauhat = obj_ripw$tauhat, ripw_se = obj_ripw$se, ripw_pval = obj_ripw$pval)
print(res)

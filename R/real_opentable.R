library("ripw")
library("tidyverse")
library("survival")
library("lfe")
library("survminer")
library("panelView")

## Preprocess data
data_full <- read.csv("../data/opentable_us_misc.csv") %>%
        filter(state.name != "District of Columbia") %>%
        select(state.name, date,
               reserv.diff,death, confirmed, Population, 
               Neighbor_soe.status, Neighbor_tot,
               soe.status, BEDS, dem_to_rep_ratio) %>%
        rename(treat = soe.status,
               state = state.name,
               vote = dem_to_rep_ratio) %>%
        mutate(date = as.Date(date),
               confirmed = log(confirmed + 1),
               beds = log(BEDS)) %>%
        select(-contains("Neighbor")) %>%
        left_join(data.frame(state = state.name, region = state.region),
                  by = "state") %>%
        na.omit()

start_date <- as.Date("2020-02-29")
end_date <- as.Date("2020-03-13")
T <- as.integer(end_date - start_date) + 1
data <- data_full %>%
    mutate(time = date - start_date) %>%
    filter(date >= start_date & date <= end_date)
nstates <- length(unique(data$state))
states <- unique(data$state)

## Summary of the variables
plot <- data %>% group_by(date) %>%
    summarize(avg_reserve_diff = mean(reserv.diff),
              avg_confirmed = mean(exp(confirmed) - 1)) %>%
    gather("measure", "value", -date) %>%
    mutate(date = as.integer(date),
           date = date - as.integer(start_date) + 1) %>%
    mutate(measure = recode(measure,
                            avg_confirmed = "Average Confirmed Cases",
                            avg_reserve_diff = "Average Reservation Difference (in %)")) %>%
    ggplot(aes(x = date, y = value)) +
    geom_line() +
    facet_wrap(measure ~ ., scales = "free_y", nrow = 2) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0.2),
                       breaks = 1:14,
                       labels = format(
                           as.Date(start_date:end_date, origin = "1970-01-01"),
                           format = "%m/%d")) +
    xlab("Date") + ylab("")
ggsave(filename = "../figs/summary_stats_1.pdf", plot, width = 6.5, height = 4.5)

plot <- data %>% group_by(state) %>%
    summarize(beds = mean(beds),
              vote = mean(vote)) %>%
    gather("measure", "value", -state) %>%
    mutate(measure = recode(measure,
                            beds = "Number of Hospital Beds",
                            vote = "Vote Share of Democrats (2016 Presidential Election)")) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 10) +
    facet_wrap(measure ~ ., scales = "free_x", nrow = 2) +
    theme_bw() +
    xlab("") + ylab("Count")
ggsave(filename = "../figs/summary_stats_2.pdf", plot, width = 6.5, height = 4.5)


####  Unweighted two-way estimator
## Center all covariates
data_reg <- data %>%
    group_by(state) %>%
    mutate(confirmed = scale(confirmed, scale = F)) %>%
    ungroup %>%
    group_by(time) %>%
    mutate(confirmed = scale(confirmed, scale = F)) %>%
    ungroup %>%
    mutate(confirmed = confirmed - mean(confirmed)) %>%
    mutate(vote = vote - mean(vote),
           beds = beds - mean(beds)) %>%
    mutate(Northeast = (region == "Northeast"),
           South = (region == "South"),
           West = (region == "West"),
           NorthCentral = (region == "North Central")) %>%
    mutate(Northeast = Northeast - mean(Northeast),
           South = South - mean(South),
           West = West - mean(West),
           NorthCentral = NorthCentral - mean(NorthCentral))
obj_unw <- felm(reserv.diff ~ treat + confirmed + treat:confirmed + treat:vote + treat:beds + treat:NorthCentral + treat:South + treat:West | factor(state) + factor(time) | 0 | 0,
                data = data_reg,
                exactDOF = TRUE,
                cmethod = "reghdfe")
unw_tauhat <- obj_unw$coefficients[1]
unw_se <- summary(obj_unw)$coefficients[1,2]
unw_pval <- 2 * pnorm(abs(unw_tauhat / unw_se), lower.tail = FALSE)
res_unw <- data.frame(
    unw_tauhat = unw_tauhat, unw_se = unw_se, unw_pval = unw_pval)
print(res_unw)

obj_unw0 <- felm(reserv.diff ~ treat + confirmed + treat:confirmed + treat:vote + treat:beds | factor(state) + factor(time) | 0 | 0,
                data = data_reg,
                exactDOF = TRUE,
                cmethod = "reghdfe")

## Test the joint significance of main terms
summary(obj_unw)
summary(obj_unw0)

## Plot the treatment paths
state_name_abb <- data.frame(state = state.name, abb = state.abb)
plot <- data %>%
    mutate(date = as.integer(date) - as.integer(start_date) + 1) %>%
    left_join(state_name_abb, by = "state") %>%
    panelview(1 ~ treat, index = c("abb", "date")) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = 1:14,
                       labels = format(
                           as.Date(start_date:end_date, origin = "1970-01-01"),
                           format = "%m/%d")) +
    xlab("Date") + ylab("State") +
    theme(legend.position = "none")
ggsave(filename = "../figs/opentable_treatment.pdf", plot, width = 8, height = 5.5)

## Show the cox-model results w/o cross-fitting
survdata <- list()
for (i in 1:length(states)){
    tmpdata <- data %>%
        filter(state == states[i]) %>%
        select(treat, time, confirmed, beds, vote, region)
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
    select(confirmed, beds, vote, region, time, time1, time2, status, id)
fit <- coxph(Surv(time1, time2, status) ~ confirmed + vote + beds + region, data = survdata)
fit0 <- coxph(Surv(time1, time2, status) ~ confirmed + vote + beds, data = survdata)
summary(fit)
summary(fit0)

## test for proportional hazard assumption
test <- cox.zph(fit)
pdf("../figs/opentable_schoenfeld.pdf", width = 9, height = 9)
ggcoxzph(test)
dev.off()

## Generate the outcomes and treatment assignments
Y <- data %>% select(state, time, reserv.diff) %>%
    spread(time, reserv.diff) %>%
    select(-state) %>%
    as.matrix
tr <- data %>% select(state, time, treat) %>%
    spread(time, treat) %>%
    select(-state) %>%
    as.matrix


## Construct covariate matrices
X_vote <- data$vote %>%
    matrix(nrow = T) %>% t %>% .[, 1]
X_beds <- data$beds %>%
    matrix(nrow = T) %>% t %>% .[, 1]
X_region <- data$region %>%
    matrix(nrow = T) %>% t %>% .[, 1] %>%
    data.frame(region = .) %>%
    model.matrix(~0+region, data = .) %>%
    .[, -1] %>% as.matrix
X_ti <- cbind(vote = X_vote, beds = X_beds, X_region)
X_confirmed <- data$confirmed %>%
    matrix(nrow = T) %>% t
X_tv <- lapply(1:T, function(t){
    tmp <- as.matrix(X_confirmed[, t])
    colnames(tmp) <- "confirmed"
    tmp
})

## Settings
pihat_other_vars <- colnames(X_ti)
muhat_other_vars <- colnames(X_ti)

## Derandomized cross-fitting
set.seed(2023)
nfolds <- 10
nreps <- 10000

pihat_params <- list(X_tv = X_tv,
                     X_ti = X_ti,
                     tic = c("confirmed", pihat_other_vars))
muhat_params <- list(X_tv = X_tv,
                     X_ti = X_ti,
                     main_tic = "confirmed",
                     int_tic = c("confirmed", pihat_other_vars),
                     joint = TRUE)
res <- ripw_twfe(Y, tr, 
                 pihat_fun = pihat_cox,
                 muhat_fun = muhat_twfe,
                 pihat_params = pihat_params,
                 muhat_params = muhat_params,
                 nreps = nreps,
                 nfolds = nfolds)
print(res)


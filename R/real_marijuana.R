library("tidyverse")
library("latex2exp")
library("survival")
library("lfe")
library("survminer")
library("stargazer")
source("ripw.R")

## Preprocess data
data_full <- read.csv("../data/pnas_processed.csv") %>%
    select(state, Year,
           ln_age_mort_rate2,
           unemployment2,
           rxdmp_original,
           rxid_original,
           pmlaw_original,
           dem_to_rep_ratio,
           TREAT) %>%
    rename(year = Year,
           mortality = ln_age_mort_rate2,
           unemployment = unemployment2,
           vote = dem_to_rep_ratio,
           treat = TREAT) %>%
    select(-contains("Neighbor")) %>%
    mutate(time = year - 1999) %>%
    na.omit()

res <- list()
yr_list <- 2002:2017

## cross-fitting
set.seed(2021)
nfolds <- 10
nstates <- length(unique(data_full$state))
foldid <- gen_cf_inds(nstates, nfolds)

for (yr in 2008:2017){
    data <- data_full %>% filter(year <= yr)
    states <- unique(data$state)

    ## Generate data for coxph
    survdata <- list()
    for (i in 1:length(states)){
        tmpdata <- data %>%
            filter(state == states[i]) %>%
            select(treat, time,
                   unemployment, vote,
                   rxdmp_original,
                   rxid_original,
                   pmlaw_original,
                   year)
        inds <- which(tmpdata$treat > 0)
        if (length(inds) > 0){
            ind <- min(inds)
            tmpdata <- tmpdata[1:ind, ] %>%
                mutate(time1 = time - 1,
                       time2 = time) %>%
                select(-time) %>%
                mutate(status = 0,
                       id = i)
            tmpdata$status[length(tmpdata$status)] <- 1
        } else {
            tmpdata <- tmpdata %>%
                mutate(time1 = time - 1,
                       time2 = time) %>%
                select(-time) %>%
                mutate(status = 0,
                       id = i)
        }
        survdata[[i]] <- tmpdata
    }
    survdata <- do.call(rbind, survdata) %>%
        select(unemployment, vote, rxdmp_original, rxid_original, pmlaw_original, year, time1, time2, status, id)

    ## test for proportional hazard assumption
    fit <- coxph(Surv(time1, time2, status) ~ unemployment + rxdmp_original + rxid_original + pmlaw_original, data = survdata)
    test <- cox.zph(fit)
    schoenfeld_pval <- as.numeric(tail(test$table[, 3], 1))

    ## Estimate generalized propensity scores with cross-fitting
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
        fit <- coxph(Surv(time1, time2, status) ~ unemployment + rxdmp_original + rxid_original + pmlaw_original, data = train)
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
            if (status == 1){
                if (length(survprob) == 1){
                    gps_est[index] <- 1 - survprob
                } else {
                    gps_est[index] <- -diff(tail(survprob, 2))
                }
            } else {
                gps_est[index] <- tail(survprob, 1)
            }
            if (status){
                adopt_time[index] <- max(survtime)
            }
        }
    }

    ## Reshaped distribution
    T <- max(data$year) - min(data$year) + 1
    adopt_time[is.na(adopt_time)] <- T + 1
    adopt_time <- T + 1 - adopt_time

    Pi <- rep(0, T + 1)
    Pi[2:T] <- 1 / 2 / T
    Pi[c(1, T+1)] <- (T + 1) / 4 / T

    EPi <- head(cumsum(rev(Pi)), -1)
    Theta <- Pi[adopt_time + 1] / gps_est

    ## Generate the outcomes and treatment assignments
    Y <- data %>% select(state, year, mortality) %>%
        spread(year, mortality) %>%
        select(-state) %>%
        as.matrix
    tr <- data %>% select(state, year, treat) %>%
        spread(year, treat) %>%
        select(-state) %>%
        as.matrix

    ## Unweighted two-way estimator     
    obj_unw <- felm(mortality ~ treat + unemployment + rxdmp_original + rxid_original + pmlaw_original | factor(state) + factor(year) | 0 | 0,
                    data = data,
                    exactDOF = TRUE,
                    cmethod = "reghdfe")
    unw_tauhat <- obj_unw$coefficients[1]
    unw_se <- summary(obj_unw)$coefficients[1,2]
    unw_pval <- 2 * pnorm(abs(unw_tauhat / unw_se), lower.tail = FALSE)

    ## RIPW estimator
    X <- data %>%
        select(treat, unemployment,
               rxdmp_original, rxid_original, pmlaw_original) %>%
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
    tmp <- yr - min(yr_list) + 1
    res[[tmp]] <- data.frame(
        year = yr,
        unw_tauhat = unw_tauhat, unw_se = unw_se, unw_pval = unw_pval,
        ripw_tauhat = obj_ripw$tauhat, ripw_se = obj_ripw$se, ripw_pval = obj_ripw$pval,
        schoenfeld_pval = schoenfeld_pval)
    print(yr)
}

res <- Reduce(rbind, res)

plot <- res %>% ggplot(aes(x = year, y = schoenfeld_pval)) +
    geom_line(size = 1.5) +
    geom_hline(yintercept = 0.05, col = "red") + 
    xlab("Year") + ylab("P-value") +
    scale_x_continuous(breaks = 2008:2017) +
    scale_y_continuous(breaks = c(0.05, 0.2, 0.4, 0.6)) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 20),
          axis.text = element_text(size = 10))
ggsave(file = "../figs/schoenfeld_marijuana.pdf", plot, width = 6, height = 5)

res_unw <- res %>%
    select(year, contains("unw")) %>%
    mutate(lb = unw_tauhat - 1.96 * unw_se,
           ub = unw_tauhat + 1.96 * unw_se) %>%
    rename(tauhat = unw_tauhat) %>%
    mutate(estimator = "UNW") %>%
    select(estimator, year, tauhat, lb, ub)
res_ripw <- res %>%
    select(year, contains("ripw")) %>%
    mutate(lb = ripw_tauhat - 1.96 * ripw_se,
           ub = ripw_tauhat + 1.96 * ripw_se) %>%
    rename(tauhat = ripw_tauhat) %>%
    mutate(estimator = "RIPW") %>%
    select(estimator, year, tauhat, lb, ub)    
plot <- rbind(res_unw, res_ripw) %>%
    ggplot(aes(x = year, y = tauhat, color = estimator, linetype = estimator)) +
    geom_line(size = 1.5) +
    geom_hline(yintercept = 0, col = "black") + 
    geom_ribbon(aes(ymin = lb, ymax = ub, color = estimator, fill = estimator, linetype = estimator), alpha = 0.25) +
    scale_x_continuous(breaks = 2008:2017,
                       expand = c(0, 0, 0, 0)) +
    scale_color_discrete(name = "Estimator") +
    scale_linetype_discrete(name = "Estimator") +
    scale_fill_discrete(name = "Estimator") +
    xlab("Year") + ylab("Estimate") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 20),
          axis.text = element_text(size = 10),
          legend.key.width = unit(1.25, "cm"))
ggsave(file = "../figs/tauhat_marijuana.pdf", plot, width = 8, height = 5)

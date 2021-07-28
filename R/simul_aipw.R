source("ripw.R")
source("utils.R")

#### Simulation setup
set.seed(20210701)
nreps <- 100
n <- 10000
T <- 4

## Parameters
sigma <- 1 # Var(\eps_{it})
sigma_X <- 1 # Var(m_{it})
sigma_tau <- 0 # Var(\tau_{it})
disX <- c(0.7, 0.3) # (P(X = 1), P(X = 2))
Xcoef <- c(0, 1, 2, 3) # the effect of covariates on Y_{it}(0) is X_i Xcoef_t
nfolds <- 10 # number of folds for cross-fitting
trprob1 <- c(0.8, 0.05, 0.05, 0.05, 0.05) # treatment probabilities when X = 1
trprob2 <- c(0.1, 0.1, 0.2, 0.3, 0.3) # treatment probabilities when X = 2

## Generate tau
tau <- gen_tau(n, T, sigma_tau = sigma_tau)
tau <- tau - mean(tau)

## Generate X
X <- sample(1:2, size = n, replace = TRUE, prob = disX)
X <- sort(X)
n1 <- sum(X == 1)
n2 <- sum(X == 2)

## Generate U
U <- sample(1:10, size = n, replace = TRUE)
U <- sort(U)

## Generate unit fixed effects
unitfe <- U * 0.5

## Generate E[Y(0)]
EY0 <- gen_EY0(n, T) + unitfe
EY0X <- center(X %*% t(Xcoef) * sigma_X)
EY0 <- EY0 + EY0X

#### Compute UNW, RIPW, and three versions of AIPW estimators
res <- list()
for (i in 1:nreps) {
    eps <- matrix(rnorm(n * T), n, T) * sigma
    treat1 <- gen_treat(n1, T, prob = trprob1)
    treat2 <- gen_treat(n2, T, prob = trprob2)
    treat <- rbind(treat1, treat2)
    Y <- EY0 + treat * tau + eps
    int_treat <- rowSums(treat) + 1    

    ## Generate folds 
    foldid <- gen_cf_inds(n, nfolds)
    
    ## Correct treatment model for UNW and RIPW
    gps_correct <- fit_gps(treat, X, foldid = foldid)

    ## Wrong treatment model for UNW and RIPW
    gps_wrong <- fit_gps(treat, rep(1, n), foldid = foldid)

    ## Correct outcome model for UNW and RIPW
    mhat_correct <- twfe_mhat(Y, treat, X,
                              interact = TRUE,
                              joint = TRUE,
                              foldid = foldid)$mhat0

    ## Wrong outcome model for UNW and RIPW
    mhat_wrong <- NULL

    ## Correct treatment model for AIPW
    mps_correct <- fit_mps(treat, X, foldid = foldid)

    ## Wrong treatment model for AIPW
    mps_wrong <- fit_mps(treat, rep(1, n), foldid = foldid)
    
    ## Correct outcome model for AIPW with fitted fe
    mhat_AIPW_correct_fe <- twfe_mhat(Y, treat, X,
                                      interact = TRUE,
                                      joint = TRUE,
                                      nfolds = 1,
                                      usefe = TRUE)

    ## Correct outcome model for AIPW without fitted fe
    mhat_AIPW_correct_nfe <- twfe_mhat(Y, treat, X,
                                       interact = TRUE,
                                       joint = TRUE,
                                       nfolds = 1,
                                       usefe = FALSE)

    ## Correct outcome model for AIPW using Mundlak regression
    mhat_AIPW_correct_mundlak <- mundlak(Y, treat, rep(X, T))
    
    ## Wrong outcome model for AIPW with fitted fe
    mhat_AIPW_wrong_fe <- twfe_mhat(Y, treat, X = NULL,
                                    interact = TRUE,
                                    joint = TRUE,
                                    nfolds = 1,
                                    usefe = TRUE)

    ## Wrong outcome model for AIPW without fitted fe
    mhat_AIPW_wrong_nfe <- twfe_mhat(Y, treat, X = NULL,
                                     interact = TRUE,
                                     joint = TRUE,
                                     nfolds = 1,
                                     usefe = FALSE)

    ## Wrong outcome model for AIPW using Mundlak regression
    mhat_AIPW_wrong_mundlak <- mundlak(Y, treat)
    
    ## Reshaped functions
    Pi_IPW <- rep(0.2, 5)
    Pi_RIPW <- c(5, 2, 2, 2, 5) / 16

    exprs <- expand.grid(out = c(TRUE, FALSE),
                         treat = c(TRUE, FALSE))
    results <- list()
    for (k in 1:nrow(exprs)){
        out_model <- exprs$out[k]
        treat_model <- exprs$treat[k]
        if (out_model){
            mhat <- mhat_correct
            mhat_AIPW_fe <- mhat_AIPW_correct_fe
            mhat_AIPW_nfe <- mhat_AIPW_correct_nfe
            mhat_AIPW_mundlak <- mhat_AIPW_correct_mundlak
        } else {
            mhat <- mhat_wrong
            mhat_AIPW_fe <- mhat_AIPW_wrong_fe
            mhat_AIPW_nfe <- mhat_AIPW_wrong_nfe
            mhat_AIPW_mundlak <- mhat_AIPW_wrong_mundlak
        }
        if (treat_model){
            gps <- gps_correct
            mps <- mps_correct
        } else {
            gps <- gps_wrong
            mps <- mps_wrong
        }
        
        Theta_UNW <- rep(1, n)
        Theta_IPW <- Pi_IPW[int_treat] / gps
        Theta_RIPW <- Pi_RIPW[int_treat] / gps
        
        est_IPW <- ripw(Y, treat, mhat, Theta_IPW)
        est_RIPW <- ripw(Y, treat, mhat, Theta_RIPW)
        est_UNW <- ripw(Y, treat, mhat, Theta_UNW)
        est_AIPW_fe <- aipw(Y, treat, mhat_AIPW_fe, mps)
        est_AIPW_nfe <- aipw(Y, treat, mhat_AIPW_nfe, mps)
        est_AIPW_mundlak <- aipw(Y, treat, mhat_AIPW_mundlak, mps)
        results[[k]] <- data.frame(
            method = c("IPW", "RIPW", "UNW", "AIPW_fe", "AIPW_nfe", "AIPW_mundlak"),
            est = c(est_IPW$tauhat, est_RIPW$tauhat,
                    est_UNW$tauhat, est_AIPW_fe$tauhat,
                    est_AIPW_nfe$tauhat,
                    est_AIPW_mundlak$tauhat),
            pval = c(est_IPW$pval, est_RIPW$pval,
                     est_UNW$pval, est_AIPW_fe$pval,
                     est_AIPW_nfe$pval,
                     est_AIPW_mundlak$pval),
            se = c(est_IPW$se, est_RIPW$se,
                   est_UNW$se, est_AIPW_fe$se,
                   est_AIPW_nfe$se,
                   est_AIPW_mundlak$se),
            out = out_model,
            treat = treat_model)
    }
    res[[i]] <- do.call(rbind, results)
    if (i %% 10 == 0) {
        print(i)
    }
}

results <- do.call(rbind, res)
save(file = "../data/simul_aipw.RData", results)

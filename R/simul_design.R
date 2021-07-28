source("ripw.R")
source("utils.R")

# Simulation setup
set.seed(20210701)
nreps <- 10000
n <- 10000
T <- 4

## Parameters
sigma <- 1 # Var(\eps_{it})
disX <- c(0.7, 0.3) # (P(X = 1), P(X = 2))
Xcoef <- c(0, 1, 2, 3) # the effect of covariates on Y_{it}(0) is X_i Xcoef_t
trprob1 <- c(0.8, 0.05, 0.05, 0.05, 0.05) # treatment probabilities when X = 1
trprob2 <- c(0.1, 0.1, 0.2, 0.3, 0.3) # treatment probabilities when X = 2

exprs <- expand.grid(vtype = c("CTE", "PTA"),
                     atype = c("constant", "uniform"))
exprs <- exprs[1:3, ]

res <- list()
for (k in 1:nrow(exprs)){
    atype <- exprs$atype[k] # type of violation
    vtype <- exprs$vtype[k] # type of a_i's
    if (vtype == "CTE"){
        sigma_X <- 0 # Var(m_{it})
        sigma_tau <- 1 # Var(\tau_{it})
    } else if (vtype == "PTA"){
        sigma_X <- 1 # Var(m_{it})
        sigma_tau <- 0 # Var(\tau_{it})
    }
    results <- data.frame(methods = NULL, est = NULL)

    ## Generate tau    
    tau <- gen_tau(n, T, sigma_tau = sigma_tau, atype = atype)
    tau <- tau - mean(tau)

    ## Generate X    
    X <- sample(1:2, size = n, replace = TRUE, prob = disX)
    X <- sort(X)
    n1 <- sum(X == 1)
    n2 <- sum(X == 2)
    X <- do.call(cbind, lapply(1:T, function(j){X}))

    ## Generate U
    U <- sample(1:10, size = n, replace = TRUE)
    U <- sort(U)

    ## Generate unit fixed effects    
    unitfe <- U * 0.5

    ## Generate E[Y(0)]
    EY0 <- gen_EY0(n, T, 0, 1) + unitfe
    EY0X <- center(t(t(X) * Xcoef * sigma_X))
    EY0 <- EY0 + EY0X

    #### Compute UNW, IPW, and RIPW estimators
    for (i in 1:nreps) {
        eps <- matrix(rnorm(n * T), n, T) * sigma
        treat1 <- gen_treat(n1, T, prob = trprob1)
        treat2 <- gen_treat(n2, T, prob = trprob2)
        treat <- rbind(treat1, treat2)
        Y <- EY0 + treat * tau + eps

        int_treat <- rowSums(treat) + 1
        gps <- c(trprob1[int_treat[1:n1]], trprob2[int_treat[n1+1:n2]])
        Pi_IPW <- rep(0.2, 5)
        Theta_IPW <- Pi_IPW[int_treat] / gps
        Pi_RIPW <- c(5, 2, 2, 2, 5) / 16
        Theta_RIPW <- Pi_RIPW[int_treat] / gps
        
        df <- data.frame(method = c("IPW", "RIPW", "UNW"))
        est_IPW <- ripw(Y, treat, NULL, Theta_IPW)
        est_RIPW <- ripw(Y, treat, NULL, Theta_RIPW)
        est_UNW <- ripw(Y, treat, NULL, Theta = rep(1, nrow(Y)))

        df$est <- c(est_IPW$tauhat, est_RIPW$tauhat, est_UNW$tauhat)
        df$pval <- c(est_IPW$pval, est_RIPW$pval, est_UNW$pval)
        results <- rbind(results, df)
        if (i %% 100 == 0) {
            print(i)
        }
    }

    results$vtype <- vtype
    results$atype <- atype
    res[[k]] <- results
}

results <- do.call(rbind, res)
save(file = "../data/simul_design.RData", results)

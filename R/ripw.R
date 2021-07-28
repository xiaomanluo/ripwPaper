source("mhat_twfe.R")

#' RIPW estimator and confidence interval
#' Inputs:
#'     Y: nXT matrix of observed outcomes
#'     tr: nXT matrix of treatment assignments
#'     muhat: nXT matrix of \hat{\mu}_{it}
#'     Theta: RIPW weights
#' Output:
#'     tauhat: RIPW estimate
#'     se: standard error
#'     tstat: t statistic
#'     pval: p-value
ripw <- function(Y, tr,
                 muhat = NULL,
                 Theta = rep(1, nrow(Y))){
    n <- nrow(Y)
    Ytd <- Y
    if (!is.null(muhat)){
        Ytd <- Ytd - muhat
    }
    Theta <- Theta / mean(Theta)
    Ytd_c <- Ytd - rowMeans(Ytd)
    tr_c <- tr - rowMeans(tr)
    
    Gamma_w <- colMeans(Theta * tr_c)
    Gamma_y <- colMeans(Theta * Ytd_c)
    Gamma_ww <- mean(Theta * rowSums(tr_c^2))
    Gamma_wy <- mean(Theta * rowSums(tr_c * Ytd_c))
    denom <- Gamma_ww - sum(Gamma_w^2)
    numer <- Gamma_wy - sum(Gamma_w * Gamma_y)
    tauhat <- numer / denom

    Ytd_c <- Ytd_c - tr_c * tauhat
    V <- (Gamma_wy + rowSums(tr_c * Ytd_c) - Ytd_c %*% Gamma_w - tr_c %*% Gamma_y) * Theta / denom
    tau_se <- sqrt(sum(V^2)) / n

    tstat <- tauhat / tau_se
    pval <- pnorm(abs(tstat), lower.tail = FALSE) * 2
    return(list(tauhat = tauhat, se = tau_se,
                tstat = tstat, pval = pval))
}

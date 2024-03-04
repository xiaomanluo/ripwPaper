##################################################################
#######
####### Helpers for similation studies in
####### Section 2.7, 5.1, & Appendix D
#######
##################################################################


#' Generate staggered assignments for simulations (Section 5.1 and Appendix D)
gen_treat <- function(n, T, prob){
    if (is.vector(prob)){
        tr <- sample(0:T, n, replace = TRUE, prob = prob)
        tr <- sapply(tr, function(x){
            vec <- rep(0, T)
            if (x == 0) {
                return(vec)
            } else {
                vec[(T + 1 - x): T] <- 1
                return(vec)
            }
        })
        tr <- t(tr)
    } else if (is.list(prob) && length(prob) == n){
        tr <- sapply(1:n, function(i){
            tmp <- sample(0:T, n, replace = TRUE, prob = prob[[i]])
            as.integer(intToBits(tmp))[T:1]
        })
        tr <- t(tr)
    }
    return(tr)
}

#' Generate E[Y(0)] for simulations (Section 5.1 and Appendix D)
gen_EY0 <- function(n, T, 
                    sigma_unitfe = 1,
                    sigma_timefe = 1,
                    sigma_m = 0,
                    nfac = 3,
                    X = NULL){
    unitfe <- center(rnorm(n)) * sigma_unitfe
    timefe <- center(rnorm(T)) * sigma_timefe
    fac1 <- matrix(runif(n * nfac), n)
    fac2 <- matrix(runif(T * nfac), T)
    if (!is.null(X)){
        fac1 <- fac1 * (X == 1)
        fac2 <- fac2 * (X == 1)
    } 
    Y0 <- fac1 %*% t(fac2) * sigma_m
    Y0 <- t(t(Y0 + unitfe) + timefe)
    return(Y0)
}

#' Generate individual treatment effects for simulations (Section 5.1 and Appendix D)
gen_tau <- function(n, T, sigma_tau = 0, atype = "uniform"){
    if (atype == "constant") {
        a <- rep(1, n)
    } else if (atype == "uniform") {
        a <- runif(n)
    } else if (atype == "gaussian") {
        a <- rnorm(n)
    }
    b <- t(center(rnorm(T)))
    tau_em <- a %*% b * sigma_tau
    tau <- matrix(tau_em, n, T)
    return(tau)
}

#' Cross fitting for generalized propensity scores for the simulation in Appendix D
fit_gps_onefold <- function(tr, X,
                            testtr = tr, Xtest = X){
    xvals <- unique(X)
    int_treat <- rowSums(tr) + 1
    gps_est <- lapply(xvals, function(x){
        table(int_treat[X == x]) / sum(X == x)
    })
    test_int_treat <- rowSums(testtr) + 1    
    n <- length(Xtest)
    gps <- rep(NA, n)
    for (i in 1:length(xvals)){
        id <- Xtest == xvals[i]
        gps[id] <- gps_est[[i]][test_int_treat[id]]
    }
    return(gps)
}

fit_gps <- function(tr, X, cf = TRUE, nfolds = 10,
                    foldid = NULL){
    if (!cf || nfolds == 1){
        gps <- fit_gps_onefold(tr, X)
        return(gps)
    } else{
        n <- nrow(tr)
        if (is.null(foldid)){
            foldid <- gen_cf_inds(n, nfolds)
        }
        gps <- rep(NA, n)
        for (i in 1:nfolds){
            inds <- foldid[[i]]
            gps[inds] <- fit_gps_onefold(tr[-inds, , drop=FALSE], X[-inds], tr[inds, , drop = FALSE], X[inds])
        }
        return(gps)
    }
}

#' Cross fitting for marginal propensity scores for the simulation in Appendix D
fit_mps <- function(tr, X, cf = TRUE, nfolds = 10,
                    foldid = NULL){
    if (is.null(foldid)){
        foldid <- gen_cf_inds(n, nfolds)
    }
    sapply(1:ncol(tr), function(j){
        fit_gps(tr[, j, drop=FALSE], X, cf, nfolds, foldid)
    })
}

#' Cross-fitting for simplified muhat estimators for the simulation in Appendix D
fit_muhat_onefold <- function(Y, tr, X = NULL,
                              Xtest = NULL,
                              ntest = NULL,
                              usefe = TRUE){
    n <- nrow(Y)
    T <- ncol(Y)
    Y <- as.numeric(Y)
    tr <- as.numeric(tr)
    
    data <- data.frame(y = Y, tr = tr)    
    data$unit <- as.factor(rep(1:n, T))
    data$time <- as.factor(rep(1:T, each = n))

    if (is.null(X)){
        formula <- "y~tr|unit+time"
    } else {
        X <- do.call(Matrix::bdiag, lapply(1:T, function(t){X}))
        X <- X[, -1]
        if (!is.null(Xtest)){
            ntest <- nrow(Xtest)
            Xtest <- do.call(Matrix::bdiag, lapply(1:T, function(t){Xtest}))
            Xtest <- Xtest[, -1]
        }
        
        xname <- paste0("x.T", 2:T)
        X <- as.data.frame(as.matrix(X))
        names(X) <- xname
        data <- data.frame(data, X)    

        formula <- paste(xname, collapse = "+")
        formula <- paste0("y~", formula)
        formula <- paste0(formula, "+tr|unit+time")
    }
    
    fit <- lfe::felm(as.formula(formula), data = data)
    coef <- coef(fit)
    tauhat <- coef["tr"]
    fe <- lfe::getfe(fit)
    unitfe <- fe$effect[1:n]
    timefe <- fe$effect[(1:T) + n]
    timefe <- timefe + mean(unitfe)
    unitfe <- unitfe - mean(unitfe)

    if (!is.null(X)){
        coef_X <- coef[xname]    
        if (is.null(Xtest)){
            muhat0 <- as.matrix(X) %*% coef_X + rep(timefe, each = n)
            if (usefe){
                muhat0 <- muhat0 + rep(unitfe, T)
            }
        } else {
            muhat0 <- Xtest %*% coef_X + rep(timefe, each = ntest)
        }
    } else {
        if (is.null(ntest)){
            muhat0 <- rep(timefe, each = n)
            if (usefe){
                muhat0 <- muhat0 + rep(unitfe, T)
            }
        } else {
            muhat0 <- rep(timefe, each = ntest)
        }
    }
    muhat1 <- muhat0 + tauhat

    return(list(muhat0 = muhat0, muhat1 = muhat1))
}

fit_muhat <- function(Y, tr, X = NULL,
                      cf = TRUE, nfolds = 10,
                      foldid = NULL, usefe = TRUE){
    n <- nrow(Y)
    T <- ncol(Y)
    
    if (!cf || nfolds == 1){
        muhat <- fit_muhat_onefold(Y, tr, X, usefe = usefe)
        muhat$muhat0 <- matrix(muhat$muhat0, ncol = T)
        muhat$muhat1 <- matrix(muhat$muhat1, ncol = T)        
        return(muhat)
    } else{
        if (!is.null(X)){
            X <- as.matrix(X)
        }
        if (is.null(foldid)){
            foldid <- gen_cf_inds(n, nfolds)
        }
        muhat0 <- matrix(NA, n, T)
        muhat1 <- matrix(NA, n, T)
        for (i in 1:nfolds){
            inds <- foldid[[i]]
            ntest <- length(inds)
            if (is.null(X)){
                muhat <- fit_muhat_onefold(
                    Y[-inds, ,drop=FALSE],
                    tr[-inds, ,drop=FALSE],
                    ntest = ntest,
                    usefe = usefe)
            } else {
                muhat <- fit_muhat_onefold(
                    Y[-inds, ,drop=FALSE],
                    tr[-inds, ,drop=FALSE],
                    X = X[-inds, ,drop=FALSE],
                    Xtest = X[inds, ,drop=FALSE],
                    ntest = ntest,
                    usefe = usefe)
            }
            muhat0[inds, ] <- matrix(muhat$muhat0, ncol = T)
            muhat1[inds, ] <- matrix(muhat$muhat1, ncol = T)
        }
        return(list(muhat0 = muhat0, muhat1 = muhat1))
    }
}

#' AIPW estimator for the simulation in Appendix D
#' 
#' Inputs:
#'     Y: nXT matrix of observed outcomes
#'     tr: nXT matrix of treatment assignments
#'     mhat: nXT matrix of \hat{m}_{it}
#'     mps: nX1 vector of marginal propensity scores
#' Output:
#'     tauhat: RIPW estimate
#'     se: standard error
#'     tstat: t statistic
#'     pval: p-value
aipw <- function(Y, tr, muhat, mps){
    n <- nrow(Y)
    T <- ncol(Y)
    if (!is.null(muhat$muhat1)){
        muhat1 <- muhat$muhat1
    } else {
        muhat1 <- matrix(0, n, T)
    }
    if (!is.null(muhat$muhat0)){
        muhat0 <- muhat$muhat0
    } else {
        muhat0 <- matrix(0, n, T)
    }
    Y1mod <- (Y - muhat1) * tr
    Y0mod <- (Y - muhat0) * (1 - tr)
    taumat <- (Y1mod - Y0mod) / mps + muhat1 - muhat0
    influence <- rowMeans(taumat)
    tauhat <- mean(influence)
    tau_se <- sd(influence) / sqrt(n)
    tstat <- tauhat / tau_se
    pval <- pnorm(abs(tstat), lower.tail = FALSE) * 2
    
    return(list(tauhat = tauhat, se = tau_se,
                tstat = tstat, pval = pval))
}

#' Simplified RIPW estimator and confidence interval for the simulation in Appendix D
#' 
#' Inputs:
#'     Y: nXT matrix of observed outcomes
#'     tr: nXT matrix of treatment assignments
#'     mhat: nXT matrix of \hat{m}_{it}
#'     Theta: RIPW weights
#' Output:
#'     tauhat: RIPW estimate
#'     se: standard error
#'     tstat: t statistic
#'     pval: p-value
ripw_simplified <- function(Y, tr, 
                            mhat = NULL,
                            Theta = rep(1, nrow(Y))){
    n <- nrow(Y)
    Ytd <- Y
    if (!is.null(mhat)){
        Ytd <- Ytd - mhat
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
    V <- (Gamma_wy - tauhat * Gamma_ww + rowSums(tr_c * Ytd_c) - Ytd_c %*% Gamma_w - tr_c %*% (Gamma_y - Gamma_w * tauhat)) * Theta / denom
    tau_se <- sd(V) / sqrt(n)

    tstat <- tauhat / tau_se
    pval <- pnorm(abs(tstat), lower.tail = FALSE) * 2
    return(list(tauhat = tauhat, se = tau_se,
                tstat = tstat, pval = pval))
}

#' Center rows and columns of vectors and matrices
center <- function(dat){
    if (is.vector(dat) || (is.matrix(dat) && ncol(dat) == 1)){
        dat <- dat - mean(dat)
    } else if (is.matrix(dat)){
        rowavg <- rowMeans(dat)
        dat <- dat - rowavg
        colavg <- colMeans(dat)
        dat <- t(t(dat) - colavg)
    } else {
        stop("dat must be a vector or a matrix")
    }
    return(dat)
}

#' Generate indices for cross-fitting
gen_cf_inds <- function(n, nfolds){
    tmp_inds <- sample(1:n, n)
    foldid <- quantile(1:n, probs = 0:nfolds / nfolds, type = 1)
    foldid[1] <- 0
    foldid <- lapply(1:nfolds, function(i){
        sort(tmp_inds[(foldid[i]+1):foldid[i+1]])
    })
    return(foldid)
}

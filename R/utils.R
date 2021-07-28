## Generate indices for cross-fitting
gen_cf_inds <- function(n, nfolds){
    tmp_inds <- sample(1:n, n)
    foldid <- quantile(1:n, probs = 0:nfolds / nfolds, type = 1)
    foldid[1] <- 0
    foldid <- lapply(1:nfolds, function(i){
        tmp_inds[(foldid[i]+1):foldid[i+1]]
    })
    return(foldid)
}

## Center rows and columns of vectors and matrices
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

## Generate staggered assignments
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

## Generate E[Y(0)]
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

## Generate individual treatment effects
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

## Cross fitting for generalized propensity scores 
fit_gps_onefold <- function(tr, X,
                            testtr = tr, testX = X){
    xvals <- unique(X)
    int_treat <- rowSums(tr) + 1
    gps_est <- lapply(xvals, function(x){
        table(int_treat[X == x]) / sum(X == x)
    })
    test_int_treat <- rowSums(testtr) + 1    
    n <- length(testX)
    gps <- rep(NA, n)
    for (i in 1:length(xvals)){
        id <- testX == xvals[i]
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
        tmp_inds <- sample(1:n, n)
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

## Cross fitting for marginal propensity scores in simul_dr.R
fit_mps <- function(tr, X, cf = TRUE, nfolds = 10,
                    foldid = NULL){
    if (is.null(foldid)){
        foldid <- gen_cf_inds(n, nfolds)
    }
    sapply(1:ncol(tr), function(j){
        fit_gps(tr[, j, drop=FALSE], X, cf, nfolds, foldid)
    })
}

## AIPW estimator for simul_dr.R
aipw <- function(Y, tr, mhat, mps){
    n <- nrow(Y)
    T <- ncol(Y)
    if (!is.null(mhat$mhat1)){
        mhat1 <- mhat$mhat1
    } else {
        mhat1 <- matrix(0, n, T)
    }
    if (!is.null(mhat$mhat0)){
        mhat0 <- mhat$mhat0
    } else {
        mhat0 <- matrix(0, n, T)
    }
    Y1mod <- (Y - mhat1) * tr
    Y0mod <- (Y - mhat0) * (1 - tr)
    taumat <- (Y1mod - Y0mod) / mps + mhat1 - mhat0
    influence <- rowMeans(taumat)
    tauhat <- mean(influence)
    tau_se <- sd(influence) / sqrt(n)
    tstat <- tauhat / tau_se
    pval <- pnorm(abs(tstat), lower.tail = FALSE) * 2
    
    return(list(tauhat = tauhat, se = tau_se,
                tstat = tstat, pval = pval))
}

## mhat using Mundlak regression without covariates
mundlak <- function(Y, tr, X = NULL){
    Wbar <- rowMeans(tr)
    WWbarc <- tr * (Wbar - mean(Wbar))
    data <- data.frame(Y = as.numeric(Y),
                       W = as.numeric(tr),
                       Wbar = Wbar,
                       WWbarc = WWbarc)
    if (!is.null(X)){
        WXbar <- rowMeans(X * tr)
        WXc <- (X - mean(X)) * tr
        WWXbarc <- (WXbar - mean(WXbar)) * tr
        data <- cbind(data, data.frame(X = X, WXbar = WXbar, WXc = WXc, WWXbarc = WWXbarc))
    }
    fit <- lm(Y ~ ., data = data)
    data0 <- data
    data0$W <- 0
    data1 <- data
    data1$W <- 1
    mhat0 <- matrix(predict(fit, newdata = data0), nrow(Y))
    mhat1 <- matrix(predict(fit, newdata = data1), nrow(Y))
    return(list(mhat0 = mhat0, mhat1 = mhat1))
}

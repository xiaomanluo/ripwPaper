## Standardize the input for two-way fixed effect regression
twfe_std_input <- function(Y, tr, X,
                           timevarying = FALSE,
                           interact = FALSE,
                           joint = FALSE,
                           unitfe = TRUE,
                           timefe = TRUE){
    ## Check the format of Y and tr
    n <- nrow(Y)
    T <- ncol(Y)
    if (nrow(tr) != n || ncol(tr) != T){
        stop("dim(Y) does not match dim(tr)")
    }

    ## Transform Y, and tr into nT X 1 vectors
    Y <- as.numeric(Y)
    tr <- as.numeric(tr)

    ## Add columns for unit and time fixed effects
    data <- data.frame(y = Y, tr = tr)    
    data$unit <- as.factor(rep(1:n, T))
    data$time <- as.factor(rep(1:T, each = n))

    ## Generate the formula for two-way regression
    fe_formula <- switch(
        as.character(2 * unitfe + timefe),
        `0` = "",
        `1` = "|time",
        `2` = "|unit",
        `3` = "|unit+time"
    )

    ## Return if no covariate is given
    if (is.null(X)){
        if (joint){
            formula <- paste0("y~tr", fe_formula)
        } else{
            formula <- paste0("y~0", fe_formula)
        }
        return(list(data = data, formula = formula,
                    X = NULL))
    }

    ## Check the format of X
    if (is.vector(X)){
        X <- as.matrix(X)
    } else if (is.list(X)){
        if (!timevarying){
            stop("Time invariant covariates must be a matrix/data.frame")
        }
        X <- lapply(X, function(x){
            if (is.numeric(x)){
                as.matrix(x)
            } else if (is.matrix(x) || is.data.frame(x)){
                x
            } else {
                stop("All elements of X must be a vector/matrix/data.frame")
            }
        })
    } else if (!is.matrix(X) && !is.data.frame(X)){
        stop("X must be a vector/matrix/data.frame or a list of them")
    }

    ## Transform X into an nT X ? matrix
    if (timevarying){
        if (is.list(X)){
            praw <- ncol(X[[1]])
        } else {
            p <- ncol(X)
            praw <- floor(p / T)
            if (p > T * praw){
                stop("The dimension of X must be divisible by the number of periods!")
            }
            X <- lapply(1:T, function(j){
                X[, (j-1)*praw + 1:praw, drop=FALSE]
            })
        }
        if (interact){
            X <- do.call(Matrix::bdiag, X)
            xname <- sapply(1:praw, function(j){
                paste0("x", j, ".T", 1:T)
            })
        } else {
            X <- do.call(rbind, X)
            xname <- paste0("x", 1:praw)
        }
    } else {
        p <- ncol(X)
        if (interact){
            X <- Matrix::bdiag(lapply(1:T, function(j){
                as.matrix(X)
            }))
            X <- X[, -(1:p), drop=FALSE]
            xname <- sapply(1:p, function(j){
                paste0("x", j, ".T", 2:T)
            })
        } else {
            xname <- paste0("x", 1:p)
        }
    }

    X <- as.data.frame(as.matrix(X))
    names(X) <- xname
    if (!unitfe && !timefe){
        X <- data.frame(Intercept = rep(1, nrow(X)), X)
    }
    data <- data.frame(data, X)

    formula <- paste(xname, collapse = "+")
    formula <- paste0("y~", formula)
    if (joint){
        formula <- paste0(formula, "+tr")
    }
    formula <- paste0(formula, fe_formula)    
    
    return(list(data = data, formula = formula,
                X = as.matrix(X)))
}

## Estimate mhat from a two-way fixed effect regression
## without cross-fitting
twfe_mhat_onefold <- function(Y, tr, X,
                              testX = NULL,
                              timevarying = FALSE,
                              interact = FALSE,
                              joint = FALSE,
                              unitfe = TRUE,
                              timefe = TRUE,
                              usefe = FALSE){
    n <- nrow(Y)
    T <- ncol(Y)
    obj <- twfe_std_input(Y, tr, X, timevarying, interact, joint, unitfe, timefe)
    if (is.null(testX)){
        overlap <- TRUE
        testX <- obj$X
    } else {
        overlap <- FALSE
        if (usefe){
            stop("Fixed effects cannot be estimated when X is not the same as testX")
        }
    }
    if (!is.null(X)){
        n_test <- floor(nrow(testX) / T)
    }

    if (!joint){
        if (usefe){
            stop("Fixed effects cannot be estimated for separate fitting")
        }
        if (is.null(X)){
            mhat1 <- matrix(0, n, T)
            mhat0 <- matrix(0, n, T)
        } else {
            fit1 <- lfe::felm(as.formula(obj$formula), data = obj$data, subset = (tr == 1))
            fit0 <- lfe::felm(as.formula(obj$formula), data = obj$data, subset = (tr == 0))
            cf1 <- as.numeric(fit1$coef)
            cf1[is.nan(cf1)] <- 0
            cf0 <- as.numeric(fit0$coef)
            cf0[is.nan(cf0)] <- 0
            mhat1 <- matrix(testX %*% cf1, n_test, T)
            mhat0 <- matrix(testX %*% cf0, n_test, T)
        }
    } else {
        fit <- lfe::felm(as.formula(obj$formula), data = obj$data)
        cf <- as.numeric(head(fit$coef, -1))
        tau <- as.numeric(tail(fit$coef, 1))
        fe <- lfe::getfe(fit)        
        timefe <- fe$effect[(1:T) + n]
        if (!overlap || !usefe){
            if (!is.null(X)){
                mu <- mean(fitted.values(fit) - tau * as.numeric(tr) - obj$X %*% cf - rep(timefe, each = n))
                mhat0 <- testX %*% cf + mu + rep(timefe, each = n_test)
                mhat0 <- matrix(mhat0, n_test)
            } else {
                mu <- mean(fitted.values(fit) - tau * as.numeric(tr) - rep(timefe, each = n))
                mhat0 <- rep(timefe, each = n) + mu
                mhat0 <- matrix(mhat0, n)
            }
        } else {
            mhat0 <- fitted.values(fit)
            if (is.null(X)){
                mhat0 <- matrix(mhat0, n)
            } else {
                mhat0 <- matrix(mhat0, n_test)
            }
        }
        mhat1 <- mhat0 + tau
    }
    mhat <- list(mhat1 = mhat1, mhat0 = mhat0)    
    return(mhat)
}

## Estimate mhat from a two-way fixed effect regression
## with cross-fitting
twfe_mhat <- function(Y, tr, X = NULL,
                      timevarying = FALSE,
                      interact = FALSE,
                      joint = FALSE,
                      unitfe = TRUE,
                      timefe = TRUE,
                      cf = TRUE, nfolds = 10,
                      foldid = NULL,
                      usefe = FALSE){
    if (!is.null(foldid)){
        nfolds <- length(foldid)
    }
    if (!cf || nfolds == 1){
        mhat <- twfe_mhat_onefold(
            Y, tr, X, NULL,
            timevarying, interact, joint,
            unitfe, timefe, usefe)
        return(mhat)
    } else {
        if (usefe){
            stop("Fixed effects cannot be estimated with cross fitting")
        }
        n <- nrow(Y)
        T <- ncol(Y)
        mhat1 <- mhat0 <- matrix(NA, n, T)
        if (!is.matrix(X)){
            X <- as.matrix(X)
        }
        if (is.null(foldid)){
            foldid <- gen_cf_inds(n, nfolds)
        }
        for (i in 1:nfolds){
            inds <- foldid[[i]]
            testX <- twfe_std_input(
                Y[inds, , drop = FALSE],
                tr[inds, , drop = FALSE],
                X[inds, , drop = FALSE],
                timevarying, interact, joint, unitfe, timefe)$X
            ## ## Need tests
            ## if (nrow(testX) == length(inds) * T){
            ##     m <- length(inds)
            ##     testX <- lapply(1:T, function(j){
            ##         testX[(j-1) * m + 1:m, , drop=FALSE]
            ##     })
            ##     testX <- Reduce(cbind, testX)
            ## } else {
            ##     stop("Dimension mismatch")
            ## }
            mhat <- twfe_mhat_onefold(
                Y[-inds, , drop = FALSE],
                tr[-inds, , drop = FALSE],
                X[-inds, , drop = FALSE],
                testX,
                timevarying, interact, joint, unitfe, timefe,
                usefe = FALSE)
            mhat1[inds, ] <- mhat$mhat1
            mhat0[inds, ] <- mhat$mhat0            
        }
        mhat <- list(mhat1 = mhat1, mhat0 = mhat0)    
        return(mhat)
    }
}

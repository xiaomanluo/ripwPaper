#### Non-linear solver for the DATE equation (Appendix B.5)

#' Solver for the DATE equation with min_{w}Pi(w) >= Pi_low
#' Inputs:
#'     design: the matrix (\tilde{w}_{(1)}, ..., \tilde{w}_{(K)})
#'     xi: weights used for the DATE
#'     init: initial value
#'     Pi_low: lower bound for min_{w}Pi(w)
#' Outputs:
#'     sol: a solution (if found)
#'     val: the value of the objective at sol
#'     init_val: the value of the objective at init
solve_date <- function(design, xi,
                       init = rep(1, ncol(design)) / ncol(design),
                       Pi_low = 0){
    T <- nrow(design)
    K <- ncol(design)
    A <- design
    JA <- A - rep(1, T) %*% t(colMeans(A))
    bmat <- A * JA - xi %*% t(diag(t(JA) %*% A))
    dateeq <- function(p){
        Ap <- A %*% p
        LHS <- bmat %*% p
        RHS <- (diag(as.numeric(Ap)) - xi %*% t(Ap)) %*% (Ap - mean(Ap))
        LHS - RHS
    }
    init_val <- dateeq(init)
    if (max(abs(init_val)) < 1e-12){
        cat("Initial value is a solution\n")
        list(sol = obj$par, val = init_val, init_val = init_val)        
    }
    normal_fac <- max(1e-6, sum(init_val^2) / 2)
    fn <- function(p){
        sum(dateeq(p)^2) / 2 / normal_fac
    }
    gr <- function(p){
        Ap <- A %*% p
        JAp <- Ap - mean(Ap)
        comp1 <- (diag(as.numeric(Ap)) - xi %*% t(Ap)) %*% JA
        comp2 <- (diag(as.numeric(JAp)) - xi %*% t(JAp)) %*% A
        Dp <- bmat - comp1 - comp2
        t(Dp) %*% dateeq(p) / normal_fac
    }
    ui1 <- rep(1, K)
    ui2 <- rep(-1, K)
    ui3 <- diag(rep(1, K))
    ui <- rbind(ui1, ui2, ui3)
    ci <- c(1-1e-4, -1-1e-4, rep(max(Pi_low - 1e-4, 0), K))

    obj <- constrOptim(init, fn, gr, ui, ci, outer.eps = 1e-6)
    val <- dateeq(obj$par)
    list(sol = obj$par, val = val, init_val = init_val)
}

#' Solver for the DATE equation that maximizes min_{w}Pi(w)
#' Inputs:
#'     design: the matrix (\tilde{w}_{(1)}, ..., \tilde{w}_{(K)})
#'     xi: weights used for the DATE
#'     init: initial value
#'     tol: tolerance for the objective value
#' Outputs:
#'     sol: a solution (if found)
#'     val: the value of the objective at sol
#'     init_val: the value of the objective at init
solve_date_uniform <- function(design, xi,
                               init = rep(1, ncol(design)) / ncol(design),
                               tol = 1e-6){
    fn <- function(Pi_low){
        obj <- solve_date(design, xi, init, Pi_low)
        max(abs(obj$val))
    }
    if (fn(0) > tol){
        cat("No solution is found\n")
        return(NULL)
    }
    left <- 0
    right <- min(init)
    if (fn(right) < tol){
        return(solve_date(design, xi, init, right))
    }
    while (right - left > 1e-4){
        mid <- (left + right) / 2
        if (fn(mid) > tol){
            right <- mid
        } else {
            left <- mid
        }
    }
    return(solve_date(design, xi, init, left))
}

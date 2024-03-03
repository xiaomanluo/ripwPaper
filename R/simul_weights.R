################################################################
#######
####### Similation study and plots in Section 2.7
#######
################################################################

library("tidyverse")
library("hrbrthemes")
library("latex2exp")
source("expr_helpers.R")

#' Treatment effect weights for weighted twfe regression
#' Inputs:
#'     tr: a matrix of treatment assignments
#'     Theta: a vector of regression weights
#' Output:
#'     wt: a matrix of weights with wt[i,t] being the weight
#'         for Y_{it}(1) - Y_{it}(0)
wls_fe_weights <- function(tr, Theta){
    n <- nrow(tr)
    Theta <- Theta / mean(Theta)
    tr_c <- tr - rowMeans(tr)
    Gamma_w <- colMeans(Theta * tr_c)
    Gamma_ww <- mean(Theta * rowSums(tr_c^2))
    denom <- Gamma_ww - sum(Gamma_w^2)
    wt <- tr_c * tr - t(t(tr) * Gamma_w)
    wt <- wt * Theta / denom
    return(wt)
}

#### Simulation setup
set.seed(20210701)
nreps <- 1000000
n <- 100
T <- 4

## Parameters
sigma <- 1 # Var(\eps_{it})
sigma_itr <- 0 # Var(m_{it})
sigma_hte <- 1 # Var(\tau_{it})
disX <- c(0.7, 0.3) # (P(X = 1), P(X = 2))
trprob1 <- trprob2 <- c(0.1, 0.2, 0.3, 0.3, 0.1)
Pi1 <- c(5, 2, 2, 2, 5) / 16 # reshaped function

## Generate X
X <- sample(1:2, size = n, replace = TRUE, prob = disX)
X <- sort(X)
n1 <- sum(X == 1)
n2 <- sum(X == 2)

#### Compute weights for the unweighted estimator
## Weights in three realizations
wt_list <- list()
for (i in 1:3){
    treat1 <- gen_treat(n1, T, prob = trprob1)
    treat2 <- gen_treat(n2, T, prob = trprob2)
    treat <- rbind(treat1, treat2)
    Theta <- rep(1, n)
    wt <- wls_fe_weights(treat, Theta)
    wt_list[[i]] <- data.frame(sample = paste("#", i),
                               wt = as.numeric(wt) * T)
}
wt_sample_unw <- do.call(rbind, wt_list)

plot_sample_unw <- wt_sample_unw %>%
    ggplot(aes(x = wt, fill = sample, color = sample)) +
    geom_density(alpha = 0.25) +
    geom_vline(xintercept = 0, linetype = "longdash") + 
    scale_x_continuous(expand = c(0, 0, 0, 0)) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    xlab("Weights") +
    ylab("Density") + 
    xlim(c(-5, 10)) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15))
ggsave("../figs/unw_sample_weights.pdf", plot_sample_unw,
       width = 5, height = 4)

## Average weights over nreps realizations
wt_avg_unw <- matrix(rep(0, n * T), n, T)
for (i in 1:nreps){
    treat1 <- gen_treat(n1, T, prob = trprob1)
    treat2 <- gen_treat(n2, T, prob = trprob2)
    treat <- rbind(treat1, treat2)
    Theta <- rep(1, n)
    wt_avg_unw <- wt_avg_unw + wls_fe_weights(treat, Theta)
}
wt_avg_unw <- data.frame(wt = as.numeric(wt_avg_unw) / nreps * T)

plot_avg_unw <- wt_avg_unw %>%
    ggplot(aes(x = wt)) +
    geom_density(alpha = 0.25, fill = "gray") +
    geom_vline(xintercept = 0, linetype = "longdash") + 
    scale_x_continuous(expand = c(0, 0, 0, 0)) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    xlab("Weights") +
    ylab("Density") +     
    xlim(c(-1, 3)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15))
ggsave("../figs/unw_avg_weights.pdf", plot_avg_unw,
       width = 5, height = 4)

#### Compute weights for the RIPW estimator
## Weights in three realizations
wt_list <- list()
for (i in 1:3){
    treat1 <- gen_treat(n1, T, prob = trprob1)
    treat2 <- gen_treat(n2, T, prob = trprob2)
    treat <- rbind(treat1, treat2)
    int_treat <- treat[, 1] + treat[, 2] + treat[, 3] + treat[, 4]
    est_gps <- table(factor(int_treat, levels = 0:4)) / n
    est_gps <- as.numeric(est_gps)
    Theta <- Pi1[int_treat + 1] / est_gps[int_treat + 1]
    wt <- wls_fe_weights(treat, Theta)
    wt_list[[i]] <- data.frame(sample = paste("#", i),
                               wt = as.numeric(wt) * T)
}
wt_sample_ripw <- do.call(rbind, wt_list)

plot_sample_ripw <- wt_sample_ripw %>%
    ggplot(aes(x = wt, fill = sample, color = sample)) +
    geom_density(alpha = 0.25) +
    geom_vline(xintercept = 0, linetype = "longdash") + 
    scale_x_continuous(expand = c(0, 0, 0, 0)) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    xlab("Weights") +
    ylab("Density") + 
    xlim(c(-8, 10.5)) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15))
ggsave("../figs/ripw_sample_weights.pdf", plot_sample_ripw,
       width = 5, height = 4)

## Average weights over nreps realizations
wt_avg_ripw <- matrix(rep(0, n * T), n, T)
for (i in 1:nreps){
    treat1 <- gen_treat(n1, T, prob = trprob1)
    treat2 <- gen_treat(n2, T, prob = trprob2)
    treat <- rbind(treat1, treat2)
    int_treat <- treat[, 1] + treat[, 2] + treat[, 3] + treat[, 4]
    est_gps <- table(factor(int_treat, levels = 0:4)) / n
    est_gps <- as.numeric(est_gps)
    Theta <- Pi1[int_treat + 1] / est_gps[int_treat + 1]
    wt_avg_ripw <- wt_avg_ripw + wls_fe_weights(treat, Theta)
}
wt_avg_ripw <- data.frame(wt = as.numeric(wt_avg_ripw) / nreps * T)

plot_avg_ripw <- wt_avg_ripw %>%
    ggplot(aes(x = wt)) +
    geom_density(alpha = 0.25, fill = "gray") +
    geom_vline(xintercept = 0, linetype = "longdash") + 
    scale_x_continuous(expand = c(0, 0, 0, 0)) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    xlab("Weights") +
    ylab("Density") +     
    xlim(c(-1, 3)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15))
ggsave("../figs/ripw_avg_weights.pdf", plot_avg_ripw,
       width = 5, height = 4)

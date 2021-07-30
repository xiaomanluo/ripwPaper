# Paper Repository

This repository contains the code to replicate all numerical results in our paper: [Double-Robust Two-Way Fixed Effect Regression For Panel Data](https://arxiv.org/abs/2107.13737).

## Introduction
All R scripts are included in the folder `R/`. The outputs and the plots are included in the folder `data/` and the folder `figs/`, respectively. 

## R scripts
The folder `R/` contains all R scripts:

- `simul_weights.R`: generate histograms of realized and average weights for unweighted and RIPW estimators (Figure 1 and 2, Section 2.5).

- `simul_design.R`, `plot_simul_design.R`: simulations on design-based inference (Figure 3, Section 5.1).

- `simul_aipw.R`, `plot_simul_aipw.R`: simulations on aggregated AIPW estimators (Figure 5, Appendix C).

- `real_marijuana`: analysis of data on medical cannabis law  (Figure 4, Section 5.2).

- `real_opentable`: analysis of OpenTable data (Section 5.3).

- `ripw.R`: compute the RIPW estimator and the confidence interval based on it.

- `mhat_twfe.R`: compute mhat given by two-way fixed effect regression.

- `date_equation.R`: non-linear solver for the DATE equation (Appendix B.5). 

- `utils.R`: helper functions. 

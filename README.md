# Paper Repository

This repository contains the code to replicate all numerical results in our paper: [Double Robust Two-Way Fixed Effect Regression For Panel Data](https://arxiv.org/abs/). It relies on the [`ripw` package](https://github.com/lihualei71/ripw) we developed to compute the RIPW-TWFE estimator.

## Introduction
All R scripts are included in the folder `R/`. The outputs and the plots are included in the folder `data/` and the folder `figs/`, respectively. To replicate all results with "one button", run the following code in Shell. It will take approximately 2 hours. 
```
chmod 755 replicate.sh
./replicate.sh
```

## R scripts
The folder `R/` contains all R scripts:

- `simul_weights.R`: generate histograms of realized and average weights for unweighted and RIPW estimators (Figure 1 and 2, Section 2.7).

- `simul_design.R`, `plot_simul_design.R`: simulations on design-based inference (Figure 3, Section 5.1).

- `simul_aipw.R`, `plot_simul_aipw.R`: simulations on aggregated AIPW estimators (Figure 6, Appendix D).

- `real_opentable`: analysis of OpenTable data (Section 5.2).


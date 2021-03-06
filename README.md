LaRA: Latent Regression Analysis
================

R package for Bayesian estimation of latent trait distributions considering hierarchical structures and partially missing covariate data. Currently `LaRA` allows the user to fit unidimensional latent regression item response models (LRMs) and their extensions for clustered observations. LRMs apply a multivariate regression equation to model the relationship between the latent trait and additional person covariates. Thus, they combine the fields of measurement models and structural analysis. LRMs are typically employed to generate plausible values in large-scale assessments.

Features
--------

-   Binary and ordered polytomous items can both be included in the analysis.
-   Population heterogeneity is taken into account either through multigroup or random intercept specifications with `mglrm()` or `rilrm()`.
-   Sampling from the posterior distribution of parameters is enriched by sampling from the full conditional distributions of missing values in person covariates.
-   Approximations for the distributions of missing values are constructed from classification and regression trees (CART).

Installing LaRA
---------------

To install the latest development version from GitHub using the `devtools` package, run:

``` r
if(!require("devtools"))
  install.packages("devtools")
devtools::install_github("jcgaasch/LaRA")
```

Dependencies
------------

LaRA relies on some routines from other R packages, where the latest CRAN version is in use: `mvtnorm`, `ucminf` and `rpart`.

Usage
-----

Imputation of five plausible values with the multigroup dataset:

``` r
library(LaRA)

## prepare data input
data(simdata_2mglrm)
Y <- simdata_2mglrm[, grep("Y", names(simdata_2mglrm), value = TRUE)]
X <- simdata_2mglrm[, grep("X", names(simdata_2mglrm), value = TRUE)]

## estimation setup: MCMC chains of length 120 with 20 initial burn-in samples
## for testing purposes (for your applications itermcmc > 10000 needed)
results <- mglrm(Y = Y, X = X, S = simdata_2mglrm$S, itermcmc = 60, burnin = 10, thin = 2)

PVs <- t(results$MCMCdraws$Theta[sample(nrow(results$MCMCdraws$Theta), size = 5), ])
```

References
----------

Aßmann, C., Gaasch, J.-C., Pohl, S., & Carstensen, C. H. (2015). Bayesian estimation in irt models with missing values in background variables. *Psychological Test and Assessment Modeling*, *54*(4), 595-618.

LaRA: Latent Regression Analysis
================

R package for Bayesian estimation of latent trait distributions considering hierarchical structures and partially missing covariate data. Currently LaRA allows the user to fit unidimensional latent regression item response models (LRMs) and their extensions for clustered observations. LRMs apply a multivariate regression equation to model the relationship between the latent trait and additional person covariates. Thus, they combine the fields of measurement models and structural analysis.

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

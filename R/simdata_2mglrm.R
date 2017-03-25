#' Test data
#'
#' Simulated dataset following the mglrm for two groups.
#' @docType data
#' @usage data(simdata_2mglrm)
#' @format A data frame with 1000 respondents allocated equally to two groups via indicator variable S. They face a test of 20 items Y of which the first 18 are binary and the last two are ordinal with 4 categories. Three covariates X are generated from a multivariate normal distribution, where the variables each have a mean of 0, unit variance and pairwise correlations equal to 0.5. X1 is transformed into a binary variable with a split value of 0. The corresponding parameters of the population model including an intercept are set to \eqn{\gamma_1=(-0.5,0.2,0.2,-0.5)'}, \eqn{\gamma_2=(1,0.4,-0.2,-0.7)'}, \eqn{\sigma_1=0.5^2} and \eqn{\sigma_2=0.6^2}. 5 percent of the observations in Y19, Y20, X1 and X2 are deleted completely at random.
"simdata_2mglrm"

#' Multilevel data
#'
#' Simulated dataset following the rilrm for ten clusters.
#' @docType data
#' @usage data(simdata_40rilrm)
#' @format A data frame with 4000 respondents nested within 40 clusters via indicator variable
#' S. They face a test of 20 items Y of which the first 18 are binary and the last two are
#' ordinal with 4 categories. Three covariates X are generated from a multivariate normal
#' distribution, where the variables each have a mean of 0, unit variance and pairwise
#' correlations equal to 0.5. X1 is transformed into a binary variable with a split value of 0.
#' The corresponding parameters of the population model including an intercept are set to
#' \eqn{\gamma=(-0.5,0.4,0.2,0.7)'}, \eqn{\sigma^2=0.6^2} and \eqn{\upsilon^2=0.5^2}. 5 percent
#' of the observations in Y18, Y20, X1 and X2 are deleted completely at random.
"simdata_40rilrm"

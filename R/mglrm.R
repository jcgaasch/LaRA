#' Multigroup Latent Regression Item Response Model
#'
#' Estimate multigroup latent regression models for binary and ordinal item response data considering partially missing covariate data.
#' @param Y data frame containing item responses. They can be binary or ordinal items. The responses must be coded starting at \code{0} or as \code{NA}. Rows of \code{Y} correspond to persons and columns correspond to items.
#' @param Ymis character string how to treat \code{NA} in \code{Y}. The default method will omit them element-wise (unbalanced panel structure) and \code{incorrect} will treat them as incorrect answers.
#' @param X data frame containing person-level covariates. They can be numeric or factor variables and contain missing values coded as \code{NA}. Rows of \code{X} correspond to persons and columns correspond to covariates. With \code{X} set to NULL (default) an interecpt-only model will be estimated.
#' @param S integer vector of individual group membership. A multigroup model with \code{S} stratifying the sample will be estimated. The most simple latent regression item response model without grouping results from \code{S} set to \code{NULL} (default).
#' @param itermcmc number of MCMC iterations.
#' @param burnin number of burnin iterations.
#' @param thin thinning interval, i.e., retain only every \code{thin}th iteration (if argument \code{thin} is used, \code{itermcmc*thin} and \code{burnin*thin} yield the total number of MCMC and burnin iterations).
#' @param tdf degrees of freedom of multivariate-t proposal distribution for category cutoff parameters for ordinal items.
#' @param cartctrl1 minimum number of observations in any terminal CART node during covariates imputation cycles.
#' @param cartctrl2 complexity parameter. Any CART split that does not decrease the overall lack of fit by a factor of \code{control2} is not attempted during covariates imputation cycles.
#' @details \code{mglrm} uses a fully Bayesian estimation approach. Independent conjugate prior distributions are chosen to develop a Metropolis-within-Gibbs sampling algorithm based on the device of data augmentation (Tanner & Wong, 1987). The function generates a sample from the posterior distribution of a one dimensional two-parameter normal ogive IRT model (Albert, 1992) \deqn{y_{ij}^*=\alpha_j \theta_i - \beta_j + \varepsilon_{ij}} including a (multivariate) regression equation of person level predictors on the mean vector of latent abilities \deqn{\theta_i=x_i\gamma_{S_i} + e_i.} Regression parameters are allowed to vary across observed groups. Partially observed person-level covariates are imputed in each sampling iteration. Sequential CART (Burgette & Reiter, 2010) are utilized as approximations to the full conditional distributions of missing values in \code{X}.
#' @return list with elements \code{MCMCdraws} and \code{M-Hacc} containing a list of posterior samples matrices (rows correspond to iterations and columns to parameters) and, if estimated, a vector of Metropolis-Hastings acceptance rates of category cutoff parameters for ordinal items.
#' @references Albert, J. H. (1992). Bayesian estimation of normal ogive item response curves using gibbs sampling. \emph{Journal of Educational Statistics}, \emph{17}(3), 251-269.
#' @references Burgette, L. F., & Reiter, J. P. (2010). Multiple imputation for missing data via sequential regression trees. \emph{American Journal of Epidemiology}, \emph{172}(9), 1070-1076.
#' @references Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior distributions by data augmentation. \emph{Journal of the American Statistical Association}, \emph{82}(398), 528-549.
#' @examples
#' ## prepare data input
#' data(simdata_2mglrm)
#' Y <- simdata_2mglrm[, grep("Y", names(simdata_2mglrm), value = TRUE)]
#' X <- simdata_2mglrm[, grep("X", names(simdata_2mglrm), value = TRUE)]
#'
#' ## estimation setup: MCMC chains of length 20 with 4 initial burn-in samples
#' ## for testing purposes (for your applications itermcmc > 10000 needed)
#' ##
#' ## not considering grouping
#' model1 <- mglrm(Y = Y, X = X, itermcmc = 10, burnin = 2, thin = 2)
#'
#' ## considering grouping
#' model2 <- mglrm(Y = Y, X = X, S = simdata_2mglrm$S, itermcmc = 10, burnin = 2, thin = 2)
#' @importFrom stats model.matrix runif rgamma rnorm pnorm qnorm predict
#' @importFrom mvtnorm rmvnorm dmvnorm rmvt dmvt
#' @importFrom ucminf ucminf
#' @importFrom rpart rpart rpart.control
#' @export
mglrm <- function(
  Y,
  Ymis = c("ignore", "incorrect"),
  X = NULL,
  S = NULL,
  itermcmc,
  burnin,
  thin = 1,
  tdf = 10,
  cartctrl1 = 5,
  cartctrl2 = 1e-04
){

  Y <- data.matrix(Y)
  YPL1 <- Y + 1
  YPL2 <- Y + 2
  N <- nrow(Y)
  J <- ncol(Y)
  Q <- apply(Y, 2, function(x){
    length(unique(x[!is.na(x)]))
  })
  Q2 <- Q - 2
  ITEMBIN <- ifelse(Q == 2, TRUE, FALSE)
  POSITEMORD <- which(Q != 2)
  if(any(is.na(Y))){
    YOBS <- !is.na(Y)
    if(match.arg(Ymis) == "incorrect"){
      Y[!YOBS] <- 0
    }
  }
  THETA <- rnorm(N)
  if(is.null(X)){
    XDM <- matrix(rep(1, N))
  } else {
    ANYXMIS <- any(is.na(X))
    if(ANYXMIS){
      XOBS <- !is.na(X)
      XMIS <- is.na(X)
      xmisord <- names(sort(colSums(XMIS)))[sort(colSums(XMIS)) > 0]
      for(k in xmisord){
        X[XMIS[, k], k] <- sample(X[XOBS[, k], k], sum(XMIS[, k]), replace = TRUE)
      }
      if(is.null(S)){
        X2IMP <- data.frame(X, THETA)
        Xcolsomit <- ncol(X2IMP)
      } else {
        X2IMP <- data.frame(X, THETA, S = as.factor(S))
        Xcolsomit <- c(ncol(X2IMP), ncol(X2IMP) - 1)
      }
    }
    XDM <- model.matrix(~., X)
  }
  KX <- ncol(XDM)
  XX <- crossprod(XDM)
  if(is.null(S)){
    S <- rep(1, N)
  }
  G <- length(unique(S))
  Ng <- table(S)
  INDG <- matrix(nrow = G, ncol = N)
  for(g in 1:G){
    INDG[g, ] <- ifelse(S == g, TRUE, FALSE)
  }
  YLAT <- matrix(0, nrow = N, ncol = J)
  GAMMA <- matrix(0, nrow = KX, ncol = G)
  SIGMA2 <- rep(1, G)
  ALPHA <- rep(1, J)
  BETA <- rep(0, J)
  XI <- t(cbind(ALPHA, BETA))
  TAU <- lapply(Q, function(x){
    if(x == 2){
      NULL
    } else {
      rep(0, x - 2)
    }
  })
  KAPPA <- lapply(Q, function(x){
    if(x == 2){
      c(-1e+05, 0, 1e+05)
    } else {
      c(-1e+05, 0, cumsum(exp(rep(0, x - 2))), 1e+05)
    }
  })
  Gamma <- matrix(0, nrow = itermcmc, ncol = G*KX)
  Sigma2 <- matrix(0, nrow = itermcmc, ncol = G)
  Alpha <- matrix(0, nrow = itermcmc, ncol = J)
  Beta <- matrix(0, nrow = itermcmc, ncol = J)
  Kappa <- matrix(0, nrow = itermcmc, ncol = sum(Q2))
  accTau <- rep(0, J)
  names(accTau) <- colnames(Y)
  Theta <- matrix(0, nrow = itermcmc, ncol = N)
  PrecGamma0 <- solve(100*diag(KX))
  shapeSigma20 <- 1
  scaleSigma20 <- 1
  scaleSigma20inv <- 1/scaleSigma20
  shapeSigma2 <- Ng/2 + shapeSigma20
  PrecXi0 <- solve(100*diag(2))
  # MCMC
  for(ii in 1:itermcmc){
    for(iii in 1:thin){
      # (1)
      for(j in 1:J){
        mu <- ALPHA[j]*THETA - BETA[j]
        FA <- pnorm(KAPPA[[j]][YPL1[, j]] - mu)
        FB <- pnorm(KAPPA[[j]][YPL2[, j]] - mu)
        YLAT[, j] <- mu + qnorm(runif(N)*(FB - FA) + FA)
      }
      # (2)
      Xitem <- cbind(THETA, -1)
      for (j in 1:J) {
        Covitem <- solve(crossprod(Xitem[YOBS[, j], ]) + PrecXi0)
        muitem <- Covitem%*%crossprod(Xitem[YOBS[, j], ], YLAT[YOBS[, j], j])
        XI[1, j] <- 0
        while(XI[1, j] <= 0){
          XI[, j] <- rmvnorm(1, muitem, Covitem)
        }
      }
      ALPHA <- XI[1, ]*(1/prod(XI[1, ]))^(1/J)
      BETA <- XI[2, ] - sum(XI[2, ])/J
      # (3)
      for(j in POSITEMORD){
        maxTau <- ucminf(par = TAU[[j]], fn = lposttau, Yj = Y[YOBS[, j], j],
          Qj = Q2[j], alpha = ALPHA[j], beta = BETA[j], Theta = THETA[YOBS[, j]],
          hessian = 1)
        hatTau <- maxTau$par
        InvhessTau <- solve(maxTau$hessian)
        TAUC <- rmvt(1, delta = hatTau, sigma = InvhessTau, df = tdf)
        ratio <- min(1, exp(
          -lposttau(TAUC, Y[YOBS[, j], j], Q2[j], ALPHA[j], BETA[j],
            THETA[YOBS[, j]]) +
          lposttau(TAU[[j]], Y[YOBS[, j], j], Q2[j], ALPHA[j], BETA[j],
            THETA[YOBS[, j]]) -
          dmvt(TAUC, delta = hatTau, sigma = InvhessTau, df = tdf, log = TRUE) +
          dmvt(TAU[[j]], delta = hatTau, sigma = InvhessTau, df = tdf, log = TRUE)
        ))
        if(runif(1) < ratio){
          TAU[[j]] <- TAUC
          KAPPA[[j]][3:Q[j]] <- cumsum(exp(TAUC))
          if(ii > burnin){
            accTau[j] <- accTau[j] + 1
          }
        }
      }
      # (4)
      for(i in 1:N){
        vartheta <- 1/(crossprod(ALPHA[YOBS[i, ]]) + 1/SIGMA2[S[i]])
        mutheta <- vartheta*(crossprod(ALPHA[YOBS[i, ]], YLAT[i, YOBS[i,]] +
          BETA[YOBS[i, ]]) + XDM[i, ]%*%GAMMA[, S[i]]/SIGMA2[S[i]])
        THETA[i] <- rnorm(1, mutheta, vartheta)
      }
      # (5)
      for (g in 1:G) {
        Covgamma <- solve(crossprod(XDM[INDG[g, ], , drop = FALSE])/SIGMA2[g] +
          PrecGamma0)
        mugamma <- Covgamma%*%crossprod(XDM[INDG[g, ], , drop = FALSE],
          THETA[INDG[g, ]])/SIGMA2[g]
        GAMMA[, g] <- rmvnorm(1, mugamma, Covgamma)
        scaleSigma2 <- 0.5*crossprod(THETA[INDG[g, ]] -
          XDM[INDG[g, ], , drop = F]%*%GAMMA[, g]) + scaleSigma20inv
        SIGMA2[g] <- 1 / rgamma(1, shape = shapeSigma2[g], rate = scaleSigma2)
      }
      if(ANYXMIS){
        # (6)
        X2IMP$THETA <- THETA
        X <- seqcart(X2IMP, xmisord, XOBS, XMIS, cartctrl1, cartctrl2)
        X <- X[, -Xcolsomit]
        XDM <- model.matrix(~., X)
        XX <- crossprod(XDM)
      }
    }
    Gamma[ii, ] <- c(GAMMA)
    Sigma2[ii, ] <- SIGMA2
    Alpha[ii, ] <- ALPHA
    Beta[ii, ] <- BETA
    Kappa[ii, ] <- unlist(lapply(KAPPA[!ITEMBIN], function(x){
      return(x[-c(1, 2, length(x))])
    }))
    Theta[ii, ] <- THETA
  }
  out <- list("MCMCdraws" = list("Gamma" = Gamma[-(1:burnin), ],
    "Sigma2" = Sigma2[-(1:burnin), ], "Alpha" = Alpha[-(1:burnin), ],
    "Beta" = Beta[-(1:burnin), ], "Kappa" = Kappa[-(1:burnin), ],
    "Theta" = Theta[-(1:burnin), ]),
    "M-Hacc" = accTau[!ITEMBIN]/(thin*itermcmc - thin*burnin))
  return(out)

}

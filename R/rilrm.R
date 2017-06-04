#' Random Intercept Latent Regression Item Response Model
#'
#' Estimate random intercept latent regression models for binary and ordinal item response data
#' considering partially missing covariate data.
#' @param Y data frame containing item responses. They can be binary or ordinal items. The
#' responses must be coded starting at \code{0} or as \code{NA}. Rows of \code{Y} correspond to
#' persons and columns correspond to items.
#' @param Ymis character string how to treat \code{NA} in \code{Y}. The default method will
#' omit them element-wise (unbalanced panel structure) and \code{incorrect} will treat them as
#' incorrect answers.
#' @param X data frame containing person-level covariates. They can be numeric or factor
#' variables and contain missing values coded as \code{NA}. Rows of \code{X} correspond to
#' persons and columns correspond to covariates. With \code{X} set to NULL (default) an
#' interecpt-only model will be estimated.
#' @param S integer vector of individual cluster membership. A random intercept model with
#' \code{S} stratifying the sample will be estimated.
#' @param itermcmc number of MCMC iterations.
#' @param burnin number of burnin iterations.
#' @param thin thinning interval, i.e., retain only every \code{thin}th iteration (if argument
#' \code{thin} is used, \code{itermcmc*thin} and \code{burnin*thin} yield the total number of
#' MCMC and burnin iterations).
#' @param tdf degrees of freedom of multivariate-t proposal distribution for category cutoff
#' parameters for ordinal items.
#' @param cartctrl1 minimum number of observations in any terminal CART node during covariates
#' imputation cycles.
#' @param cartctrl2 complexity parameter. Any CART split that does not decrease the overall
#' lack of fit by a factor of \code{control2} is not attempted during covariates imputation
#' cycles.
#' @details \code{rilrm} uses a fully Bayesian estimation approach. Independent conjugate prior
#' distributions are chosen to develop a Metropolis-within-Gibbs sampling algorithm based on
#' the device of data augmentation (Tanner & Wong, 1987). The function generates a sample from
#' the posterior distribution of a one dimensional two-parameter normal ogive IRT model
#' (Albert, 1992) \deqn{y_{cij}^*=\alpha_j \theta_{ci} - \beta_j + \varepsilon_{cij}} including
#' a (multivariate) regression equation of a cluster level random intercept and person level
#' predictors on the mean vector of latent abilities
#' \deqn{\theta_{ci}=\omega_c + x_{ci}\gamma + e_{ci}.} This corresponds to the most basic
#' multilevel specification (Fox & Glas, 2001). Partially observed person-level covariates are
#' imputed in each sampling iteration. Sequential CART (Burgette & Reiter, 2010) are utilized
#' as approximations to the full conditional distributions of missing values in \code{X}.
#' @return list with elements \code{MCMCdraws} and \code{M-Hacc} containing a list of posterior
#' samples matrices (rows correspond to iterations and columns to parameters) and, if
#' estimated, a vector of Metropolis-Hastings acceptance rates of category cutoff parameters
#' for ordinal items.
#' @references Albert, J. H. (1992). Bayesian estimation of normal ogive item response curves
#' using gibbs sampling. \emph{Journal of Educational Statistics}, \emph{17}(3), 251-269.
#' @references Burgette, L. F., & Reiter, J. P. (2010). Multiple imputation for missing data
#' via sequential regression trees. \emph{American Journal of Epidemiology}, \emph{172}(9),
#' 1070-1076.
#' @references Fox, J.-P., & Glas, C. A. W. (2001). Bayesian estimation of a multilevel irt
#' model using Gibbs sampling. \emph{Psychometrika}, \emph{66}(2), 271-288.
#' @references Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior distributions
#' by data augmentation. \emph{Journal of the American Statistical Association},
#' \emph{82}(398), 528-549.
#' @examples
#' ## prepare data input
#' data(simdata_40rilrm)
#' Y <- simdata_40rilrm[, grep("Y", names(simdata_40rilrm), value = TRUE)]
#' X <- simdata_40rilrm[, grep("X", names(simdata_40rilrm), value = TRUE)]
#'
#' ## estimation setup: MCMC chains of length 20 with 4 initial burn-in samples
#' ## for testing purposes (for your applications itermcmc > 10000 needed)
#' ##
#' results <- rilrm(Y = Y, X = X, S = simdata_40rilrm$S, itermcmc = 10, burnin = 2, thin = 2)
#' @importFrom stats model.matrix runif rgamma rnorm pnorm qnorm predict
#' @importFrom mvtnorm rmvnorm dmvnorm rmvt dmvt
#' @importFrom ucminf ucminf
#' @importFrom rpart rpart rpart.control
#' @export
rilrm <- function(
  Y,
  Ymis = c("ignore", "incorrect"),
  X = NULL,
  S,
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
  QMI2 <- Q - 2
  ITEMBIN <- ifelse(Q == 2, TRUE, FALSE)
  POSITEMORD <- which(Q != 2)
  YOBS <- !is.na(Y)
  if(match.arg(Ymis) == "incorrect"){
    Y[!YOBS] <- 0
  }
  suS <- sort(unique(S))
  CL <- length(unique(suS))
  for(cl in 1:CL){
    S[S == suS[cl]] <- cl
  }
  Ncl <- table(S)
  Dcl <- matrix(0, nrow = N, ncol = CL)
  for(i in 1:N){
    Dcl[i, S[i]] <- 1
  }
  DD <- crossprod(Dcl)
  OMEGA <- rnorm(CL)
  OMEGAi <- Dcl%*%OMEGA
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
      X2IMP <- data.frame(X, THETA, OMEGAi)
      Xcolsomit <- c(ncol(X2IMP), ncol(X2IMP) - 1)
    }
    XDM <- model.matrix(~., X)
  }
  KX <- ncol(XDM)
  XX <- crossprod(XDM)
  YLAT <- matrix(0, nrow = N, ncol = J)
  GAMMA <- matrix(rep(0, KX), ncol = 1)
  SIGMA2 <- 1
  UPSILON2 <- 1
  ALPHA <- rep(1, J)
  BETA <- rep(0, J)
  XI <- rbind(ALPHA, BETA)
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
  Omega <- matrix(0, nrow = itermcmc, ncol = CL)
  Theta <- matrix(0, nrow = itermcmc, ncol = N)
  Gamma <- matrix(0, nrow = itermcmc, ncol = KX)
  Sigma2 <- matrix(0, nrow = itermcmc, ncol = 1)
  Upsilon2 <- matrix(0, nrow = itermcmc, ncol = 1)
  Alpha <- matrix(0, nrow = itermcmc, ncol = J)
  Beta <- matrix(0, nrow = itermcmc, ncol = J)
  Kappa <- matrix(0, nrow = itermcmc, ncol = sum(QMI2))
  accTau <- rep(0, J)
  names(accTau) <- colnames(Y)
  PrecGamma0 <- solve(100*diag(KX))
  shapeSigma20 <- 1
  scaleSigma20 <- 1
  scaleSigma20inv <- 1/scaleSigma20
  shapeSigma2 <- N/2 + shapeSigma20
  shapeUpsilon20 <- 1
  scaleUpsilon20 <- 1
  scaleUpsilon20inv <- 1/scaleUpsilon20
  shapeUpsilon2 <- CL/2 + shapeUpsilon20
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
        maxTau <- ucminf(par = TAU[[j]], fn = lposttau, Yj = Y[YOBS[, j], j], Qj = QMI2[j],
          alpha = ALPHA[j], beta = BETA[j], Theta = THETA[YOBS[, j]], hessian = 1)
        hatTau <- maxTau$par
        InvhessTau <- solve(maxTau$hessian)
        TAUC <- rmvt(1, delta = hatTau, sigma = InvhessTau, df = tdf)
        ratio <- min(1, exp(
          -lposttau(TAUC, Y[YOBS[, j], j], QMI2[j], ALPHA[j], BETA[j], THETA[YOBS[, j]]) +
          lposttau(TAU[[j]], Y[YOBS[, j], j], QMI2[j], ALPHA[j], BETA[j], THETA[YOBS[, j]]) -
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
        vartheta <- 1/(crossprod(ALPHA[YOBS[i, ]]) + 1/SIGMA2)
        mutheta <- vartheta*(crossprod(ALPHA[YOBS[i, ]], YLAT[i, YOBS[i,]] + BETA[YOBS[i, ]]) +
          (OMEGAi[i] + XDM[i, ]%*%GAMMA)/SIGMA2)
        THETA[i] <- rnorm(1, mutheta, sqrt(vartheta))
      }
      # (5)
      Covgamma <- solve(XX/SIGMA2 + PrecGamma0)
      mugamma <- Covgamma%*%crossprod(XDM, THETA - OMEGAi)/SIGMA2
      GAMMA <- t(rmvnorm(1, mugamma, Covgamma))
      XGAMMA <- XDM%*%GAMMA
      # (6)
      scaleSigma2 <- 0.5*crossprod(THETA - OMEGAi - XGAMMA) + scaleSigma20inv
      SIGMA2 <- 1/rgamma(1, shape = shapeSigma2, rate = scaleSigma2)
      # (7)
      Covri <- solve(DD/SIGMA2 + (1/UPSILON2)*diag(CL))
      muri <- Covri%*%(crossprod(Dcl, THETA - XGAMMA)/SIGMA2)
      OMEGA <- c(rmvnorm(1, muri, Covri))
      OMEGAi <- Dcl%*%OMEGA
      # (8)
      scaleUpsilon2 <- 0.5*crossprod(OMEGA) + scaleUpsilon20inv
      UPSILON2 <- 1/rgamma(1, shape = shapeUpsilon2, rate = scaleUpsilon2)
      if(ANYXMIS){
        # (9)
        X2IMP$THETA <- THETA
        X2IMP$OMEGAi <- OMEGAi
        X <- seqcart(X2IMP, xmisord, XOBS, XMIS, cartctrl1, cartctrl2)
        X <- X[, -Xcolsomit]
        XDM <- model.matrix(~., X)
        XX <- crossprod(XDM)
      }
    }
    Omega[ii, ] <- OMEGA
    Theta[ii, ] <- THETA
    Gamma[ii, ] <- GAMMA
    Sigma2[ii, ] <- SIGMA2
    Upsilon2[ii, ] <- UPSILON2
    Alpha[ii, ] <- ALPHA
    Beta[ii, ] <- BETA
    Kappa[ii, ] <- unlist(lapply(KAPPA[!ITEMBIN], function(x){
      return(x[-c(1, 2, length(x))])
    }))
  }
  bi <- 1:burnin
  out <- list("MCMCdraws" = list("Omega" = Omega[-bi, ],
    "Theta" = Theta[-bi, ], "Gamma" = Gamma[-bi, ],
    "Sigma2" = Sigma2[-bi, ], "Upsilon2" = Upsilon2[-bi, ],
    "Alpha" = Alpha[-bi, ], "Beta" = Beta[-bi, ],
    "Kappa" = Kappa[-bi, ]),
    "M-Hacc" = accTau[!ITEMBIN]/(thin*itermcmc - thin*burnin)
  )
  return(out)

}

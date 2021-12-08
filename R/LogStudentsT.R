################################################################################
################# IMPLEMENTATION OF THE LOG-STUDENT'S T MODEL ##################
################################################################################

#' @title MCMC algorithm for the log-student's t model
#' @description  Adaptive Metropolis-within-Gibbs algorithm with univariate
#'     Gaussian random walk proposals for the log-student's T model (no mixture)
#' @inheritParams MCMC_LN
#' @param Q Update period for the \eqn{\lambda_{i}}’s
#' @param ar Optimal acceptance rate for the adaptive Metropolis-Hastings
#'     updates
#' @param nu0 Starting value for \eqn{v}. If not provided, then it will be
#'     randomly generated from a gamma distribution.
#' @return A matrix with \eqn{N / thin + 1} rows. The columns are the MCMC
#'     chains for \eqn{\beta} (\eqn{k} columns), \eqn{\sigma^2} (1 column),
#'     \eqn{\theta} (1 column, if appropriate), \eqn{\lambda} (\eqn{n} columns,
#'     not provided for log-normal model), \eqn{\log(t)} (\eqn{n} columns,
#'     simulated via data augmentation) and the logarithm of the adaptive
#'     variances (the number varies among models). The latter allows the user to
#'     evaluate if the adaptive variances have been stabilized.
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations.
#'
#' LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#'
#' @export
MCMC_LST <- function(N,
                     thin,
                     burn,
                     Time,
                     Cens,
                     X,
                     Q = 1,
                     beta0 = NULL,
                     sigma20 = NULL,
                     nu0 = NULL,
                     prior = 2,
                     set = TRUE,
                     eps_l = 0.5,
                     eps_r = 0.5,
                     ar = 0.44) {

  # Sample starting values if not given
  if (is.null(beta0)) beta0 <- beta.sample(n = ncol(X))
  if (is.null(sigma20)) sigma20 <- sigma2.sample()
  if (is.null(nu0)) nu0 <- nu.sample()

  MCMC.param.check(N,
                   thin,
                   burn,
                   Time,
                   Cens,
                   X,
                   beta0,
                   sigma20,
                   prior,
                   set,
                   eps_l,
                   eps_r)



  k <- length(beta0)
  n <- length(Time)
  N.aux <- round(N / thin, 0)
  if (prior == 1) {
    p <- 1 + k / 2
  }
  if (prior == 2) {
    p <- 1
  }

  beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
  beta[1, ] <- beta0
  sigma2 <- rep(0, times = N.aux + 1)
  sigma2[1] <- sigma20
  nu <- rep(0, times = N.aux + 1)
  nu[1] <- nu0
  logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
  logt[1, ] <- log(Time)
  lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
  lambda[1, ] <- stats::rgamma(n, shape = nu0 / 2, rate = nu0 / 2)
  accept.nu <- 0
  pnu.aux <- 0
  ls.nu <- rep(0, times = N.aux + 1)

  beta.aux <- beta[1, ]
  sigma2.aux <- sigma2[1]
  nu.aux <- nu[1]
  logt.aux <- logt[1, ]
  lambda.aux <- lambda[1, ]
  ls.nu.aux <- ls.nu[1]

  i_batch <- 0

  for (iter in 2:(N + 1)) {
    i_batch <- i_batch + 1

    Lambda <- diag(lambda.aux)
    AUX1 <- (t(X) %*% Lambda %*% X)
    if (det(AUX1) != 0) {
      AUX <- solve(AUX1)
      mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
      Sigma.aux <- sigma2.aux * AUX
      beta.aux <- mvrnormArma(n = 1, mu = mu.aux, Sigma = Sigma.aux)
    }

    shape.aux <- (n + 2 * p - 2) / 2
    rate.aux <- 0.5 * t(logt.aux - X %*% beta.aux) %*% Lambda %*%
      (logt.aux - X %*% beta.aux)
    if (rate.aux > 0 & is.na(rate.aux) == FALSE) {
      sigma2.aux <- (stats::rgamma(1, shape = shape.aux,
                                   rate = rate.aux)) ^ (-1)
    }

    MH.nu <- MH_nu_LST(omega2 = exp(ls.nu.aux),
                       beta.aux,
                       lambda.aux,
                       nu.aux,
                       prior)
    nu.aux <- MH.nu$nu
    accept.nu <- accept.nu + MH.nu$ind
    pnu.aux <- pnu.aux + MH.nu$ind

    if ((iter - 1) %% Q == 0 && iter - 1 <= burn) {
      shape1.aux <- (nu.aux + 1) / 2
      rate1.aux <- 0.5 * (nu.aux + ((logt.aux - X %*% beta.aux) ^ 2) /
                            sigma2.aux)
      lambda.aux <- stats::rgamma(n,
                                  shape = rep(shape1.aux, times = n),
                                  rate = rate1.aux)
    }

    logt.aux <- logt.update.SMLN(Time,
                                 Cens,
                                 X,
                                 beta.aux,
                                 sigma2.aux / lambda.aux,
                                 set,
                                 eps_l,
                                 eps_r)

    if (i_batch == 50) {
      pnu.aux <- pnu.aux / 50
      Pnu.aux <- as.numeric(pnu.aux < ar)
      ls.nu.aux <- ls.nu.aux + ((-1) ^ Pnu.aux) *
        min(0.01, 1 / sqrt(iter))
      i_batch <- 0
      pnu.aux <- 0
    }

    if (iter %% thin == 0) {
      beta[iter / thin + 1, ] <- beta.aux
      sigma2[iter / thin + 1] <- sigma2.aux
      nu[iter / thin + 1] <- nu.aux
      logt[iter / thin + 1, ] <- logt.aux
      lambda[iter / thin + 1, ] <- lambda.aux
      ls.nu[iter / thin + 1] <- ls.nu.aux
    }
    if ((iter - 1) %% 1e+05 == 0) {
      cat(paste("Iteration :", iter, "\n"))
    }
  }

  cat(paste("AR nu :", round(accept.nu / N, 2), "\n"))

  chain <- cbind(beta, sigma2, nu, lambda, logt, ls.nu)

  beta.cols <- paste("beta.", seq(ncol(beta)), sep = "")
  lambda.cols <- paste("lambda.", seq(ncol(lambda)), sep = "")
  logt.cols <- paste("logt.", seq(ncol(logt)), sep = "")

  colnames(chain) <- c(beta.cols,
                       "sigma2",
                       "nu",
                       lambda.cols,
                       logt.cols,
                       "ls.nu")


  if (burn > 0) {
    burn.period <- 1:(burn / thin)
    chain <- chain [-burn.period, ]
  }

  return(chain)
}

#' @title Log-marginal Likelihood estimator for the log-student's t model
#' @inheritParams MCMC_LN
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @param Q Update period for the \eqn{\lambda_{i}}’s
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations.
#'
#' LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#'
#' LST.LML <- LML_LST(thin = 20, Time = cancer[, 1], Cens = cancer[, 2],
#'                    X = cancer[, 3:11], chain = LST)
#'
#' @export
LML_LST <- function(thin,
                    Time,
                    Cens,
                    X,
                    chain,
                    Q = 1,
                    prior = 2,
                    set = TRUE,
                    eps_l = 0.5,
                    eps_r = 0.5) {
  chain <- as.matrix(chain)
  n <- length(Time)
  N <- dim(chain)[1]
  k <- dim(X)[2]
  if (prior == 1) {
    p <- 1 + k / 2
  }
  if (prior == 2) {
    p <- 1
  }
  omega2.nu <- exp(stats::median(chain[, k + 3 + 2 * n]))

  chain.nonadapt <- MCMC.LST.NonAdapt(N = N * thin,
                                      thin,
                                      Q,
                                      Time,
                                      Cens,
                                      X,
                                      beta0 = t(chain[N, 1:k]),
                                      sigma20 = chain[N, k + 1],
                                      nu0 = chain[N, k + 2],
                                      prior,
                                      set,
                                      eps_l,
                                      eps_r,
                                      omega2.nu)
  chain.nonadapt <- chain.nonadapt[-1, ]
  sigma2.star <- stats::median(chain.nonadapt[, k + 1])
  nu.star <- stats::median(chain.nonadapt[, k + 2])
  if (k > 1) {
    beta.star <- apply(chain.nonadapt[, 1:k], 2, "median")
  } else {
    beta.star <- stats::median(chain.nonadapt[, 1])
  }

  # Log-LIKELIHOOD ORDINATE
  LL.ord <- log_lik_LST(Time,
                        Cens,
                        X,
                        beta = beta.star,
                        sigma2 = sigma2.star,
                        nu = nu.star,
                        set,
                        eps_l,
                        eps_r)
  cat("Likelihood ordinate ready!\n")

  # PRIOR ORDINATE
  LP.ord <- prior_LST(beta = beta.star,
                      sigma2 = sigma2.star,
                      nu = nu.star,
                      prior,
                      logs = TRUE)
  cat("Prior ordinate ready!\n")


  logt0 <- t(chain.nonadapt[N, (n + k + 3):(2 * n + k + 2)])
  lambda0 <- t(chain.nonadapt[N, (k + 3):(k + 2 + n)])

  chain.sigma2 <- MCMCR.nu.LST(N = N * thin,
                               thin = thin,
                               Q,
                               Time,
                               Cens,
                               X,
                               beta0 = t(chain.nonadapt[N, 1:k]),
                               sigma20 = chain.nonadapt[N, (k + 1)],
                               nu0 = nu.star,
                               logt0 = logt0,
                               lambda0 = lambda0,
                               prior = prior,
                               set,
                               eps_l,
                               eps_r)
  cat("Reduced chain.sigma2 ready!\n")

  # POSTERIOR ORDINATE - nu
  # USING AN ADAPTATION of the CHIB AND JELIAZKOV (2001) METHODOLOGY
  po1.nu <- rep(0, times = N)
  po2.nu <- rep(0, times = N)
  for (i in 1:N) {
    lambda.po1 <- t(chain.nonadapt[i, (k + 3):(k + 2 + n)])
    mean <- as.numeric(chain.nonadapt[i, (k + 2)])
    po1.nu[i] <- alpha_nu(nu0 = as.numeric(chain.nonadapt[i, (k + 2)]),
                          nu1 = nu.star,
                          lambda = lambda.po1,
                          prior = prior) * stats::dnorm(x = nu.star,
                                                        mean = mean,
                                                        sd = sqrt(omega2.nu))
    nu.aux <- stats::rnorm(n = 1, mean = nu.star, sd = sqrt(omega2.nu))

    lambda.po2 <- t(chain.sigma2[i + 1, (k + 2):(k + 1 + n)])
    po2.nu[i] <- alpha_nu(nu0 = nu.star,
                          nu1 = nu.aux,
                          lambda = lambda.po2,
                          prior = prior)
  }
  PO.nu <- mean(po1.nu) / mean(po2.nu)
  cat("Posterior ordinate nu ready!\n")

  # POSTERIOR ORDINATE - sigma 2
  shape <- (n + 2 * p - 2) / 2
  po.sigma2 <- rep(0, times = N)
  for (i in 1:N) {
    aux1 <- (chain.sigma2[i, (k + 2 + n):(k + 1 + 2 * n)]) - X %*%
      (chain.sigma2[i, 1:k])
    aux2 <- diag(as.vector(t(chain.sigma2[i, (k + 2):(k + 1 + n)])))
    rate.aux <- as.numeric(0.5 * t(aux1) %*% aux2 %*% aux1)
    po.sigma2[i] <- MCMCpack::dinvgamma(sigma2.star,
                                        shape = shape,
                                        scale = rate.aux)
  }
  PO.sigma2 <- mean(po.sigma2)
  cat("Posterior ordinate sigma2 ready!\n")

  # POSTERIOR ORDINATE - beta
  logt0 <- t(chain.nonadapt[N, (n + k + 3):(2 * n + k + 2)])
  lambda0 <- t(chain.nonadapt[N, (k + 3):(k + 2 + n)])

  chain.beta <- MCMCR.sigma2.nu.LST(N = N * thin,
                                    thin = thin,
                                    Q,
                                    Time,
                                    Cens,
                                    X,
                                    beta0 = t(chain.nonadapt[N, 1:k]),
                                    sigma20 = sigma2.star,
                                    nu0 = nu.star,
                                    logt0 = logt0,
                                    lambda0 = lambda0,
                                    prior,
                                    set,
                                    eps_l,
                                    eps_r)

  cat("Reduced chain.beta ready!\n")
  po.beta <- rep(0, times = N)
  for (i in 1:N) {
    aux0.beta <- diag(as.vector(t(chain.beta[(i + 1), (k + 1):(k + n)])))
    aux1.beta <- solve(t(X) %*% aux0.beta %*% X)
    aux2.beta <- aux1.beta %*% t(X) %*% aux0.beta
    mu.aux <- as.vector(aux2.beta %*%
                          ((chain.beta[(i + 1), (k + 1 + n):(k + 2 * n)])))
    po.beta[i] <- mvtnorm::dmvnorm(beta.star,
                                   mean = mu.aux,
                                   sigma = sigma2.star * aux1.beta,
                                   log = FALSE)
  }
  PO.beta <- mean(po.beta)
  cat("Posterior ordinate beta ready!\n")

  # TAKING LOGARITHM
  LPO.nu <- log(PO.nu)
  LPO.sigma2 <- log(PO.sigma2)
  LPO.beta <- log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML <- LL.ord + LP.ord - LPO.nu - LPO.sigma2 - LPO.beta

  return(list(LL.ord = LL.ord,
              LP.ord = LP.ord,
              LPO.nu = LPO.nu,
              LPO.sigma2 = LPO.sigma2,
              LPO.beta = LPO.beta,
              LML = LML))
}

#' @title Deviance information criterion for the log-student's t model
#' @description Deviance information criterion is based on the deviance function
#'     \eqn{D(\theta, y) = -2 log(f(y|\theta))} but also incorporates a
#'     penalization factor of the complexity of the model
#' @inheritParams MCMC_LN
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations.
#'
#' LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' LST.DIC <- DIC_LST(Time = cancer[, 1], Cens = cancer[, 2],
#'                    X = cancer[, 3:11], chain = LST)
#'
#' @export
DIC_LST <- function(Time,
                    Cens,
                    X,
                    chain,
                    set = TRUE,
                    eps_l = 0.5,
                    eps_r = 0.5) {
  chain <- as.matrix(chain)
  N <- dim(chain)[1]
  k <- dim(X)[2]
  LL <- rep(0, times = N)

  for (iter in 1:N) {
    LL[iter] <- log_lik_LST(Time,
                            Cens,
                            X,
                            beta = as.vector(chain[iter, 1:k]),
                            sigma2 = chain[iter, k + 1],
                            nu = chain[iter, k + 2],
                            set,
                            eps_l,
                            eps_r)
  }

  aux <- apply(chain[, 1:(k + 2)], 2, "median")
  pd <- -2 * mean(LL) + 2 * log_lik_LST(Time,
                                        Cens,
                                        X,
                                        beta = aux[1:k],
                                        sigma2 = aux[k + 1],
                                        nu = aux[k + 2],
                                        set,
                                        eps_l,
                                        eps_r)
  pd.aux <- k + 2

  DIC <- -2 * mean(LL) + pd

  cat(paste("Effective number of parameters :", round(pd, 2), "\n"))
  cat(paste("Actual number of parameters :", pd.aux, "\n"))
  return(DIC)
}

#' @title Case deletion analysis for the log-student's t model
#' @description Leave-one-out cross validation analysis. The function returns a
#'     matrix with n rows. The first column contains the logarithm of the CPO
#'     (Geisser and Eddy, 1979). Larger values of the CPO indicate better
#'     predictive accuracy of the model. The second and third columns contain
#'     the KL divergence between \eqn{\pi(\beta, \sigma^2,  \theta | t_{-i})}
#'     and \eqn{\pi(\beta, \sigma^2,  \theta | t)} and its calibration index
#'     \eqn{p_i}, respectively.
#' @inheritParams MCMC_LN
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations.
#'
#' LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' LST.CD <- CaseDeletion_LST(Time = cancer[, 1], Cens = cancer[, 2],
#'                            cancer[, 3:11], chain = LST)
#'
#' @export
CaseDeletion_LST <- function(Time,
                             Cens,
                             X,
                             chain,
                             set = TRUE,
                             eps_l = 0.5,
                             eps_r = 0.5) {
  chain <- as.matrix(chain)
  n <- dim(X)[1]
  k <- dim(X)[2]
  logCPO <- rep(0, times = n)
  KL.aux <- rep(0, times = n)
  N <- dim(chain)[1]

  for (s in 1:n) {
    aux1 <- rep(0, times = N)
    aux2 <- rep(0, times = N)
    for (ITER in 1:N) {
      aux2[ITER] <- log.lik.LST(Time[s],
                                Cens[s],
                                X[s, ],
                                beta = as.vector(chain[ITER, 1:k]),
                                sigma2 = chain[ITER, k + 1],
                                nu = chain[ITER, k + 2],
                                set,
                                eps_l,
                                eps_r)
      aux1[ITER] <- exp(-aux2[ITER])
    }
    logCPO[s] <- -log(mean(aux1))
    KL.aux[s] <- mean(aux2)
    if (KL.aux[s] <- logCPO[s] < 0) {
      cat(paste("Numerical problems for observation:", s, "\n"))
    }
  }
  KL <- abs(KL.aux - logCPO)
  CALIBRATION <- 0.5 * (1 + sqrt(1 - exp(-2 * KL)))
  return(cbind(logCPO, KL, CALIBRATION))
}

#' @title Outlier detection for observation for the log-student's t model
#' @description This returns a unique number corresponding to the Bayes Factor
#'     associated to the test \eqn{M_0: \Lambda_{obs} = \lambda_{ref}} versus
#'     \eqn{M_1: \Lambda_{obs}\neq \lambda_{ref}} (with all other
#'     \eqn{\Lambda_j,\neq obs} free). The value of \eqn{\lambda_{ref}} is
#'     required as input. The user should expect long running times for the
#'     log-Student’s t model, in which case a reduced chain given
#'     \eqn{\Lambda_{obs} = \lambda_{ref}} needs to be generated
#' @inheritParams MCMC_LN
#' @param Q Update period for the \eqn{\lambda_{i}}’s
#' @param burn Burn-in period
#' @param ref Reference value \eqn{\lambda_{ref}} or \eqn{u_{ref}}
#' @param obs Indicates the number of the observation under analysis
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @param ar Optimal acceptance rate for the adaptive Metropolis-Hastings
#'     updates
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations.
#'
#' LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#'
#' LST.Outlier <- BF_lambda_obs_LST(N = 100, thin = 20 , burn = 1, ref = 1,
#'                                  obs = 1, Time = cancer[, 1],
#'                                  Cens = cancer[, 2], X = cancer[, 3:11],
#'                                  chain = LST)
#'
#' @export
BF_lambda_obs_LST <- function(N,
                              thin,
                              burn,
                              ref,
                              obs,
                              Time,
                              Cens,
                              X,
                              chain,
                              Q = 1,
                              prior = 2,
                              set = TRUE,
                              eps_l = 0.5,
                              eps_r =0.5,
                              ar = 0.44) {
  chain <- as.matrix(chain)
  aux1 <- Post_lambda_obs_LST(obs, ref, X, chain)
  aux2 <- CFP.obs.LST(N,
                      thin,
                      Q,
                      burn,
                      ref,
                      obs,
                      Time,
                      Cens,
                      X,
                      chain,
                      prior,
                      set,
                      eps_l,
                      eps_r,
                      ar = 0.44)
  return(aux1 * aux2)
}

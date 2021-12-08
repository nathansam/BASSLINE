################################################################################
############### IMPLEMENTATION OF THE LOG-EXPONENTIAL POWER MODEL ##############
################################################################################
#' @title MCMC algorithm for the log-exponential power model
#' @description  Adaptive Metropolis-within-Gibbs algorithm with univariate
#'     Gaussian random walk proposals for the log-exponential model
#' @inheritParams MCMC_LN
#' @param alpha0 Starting value for \eqn{\alpha}. If not provided, then it will
#'     be randomly generated from a uniform distribution.
#' @param ar Optimal acceptance rate for the adaptive Metropolis-Hastings
#'     updates
#' @return A matrix with \eqn{N / thin + 1} rows. The columns are the MCMC
#'     chains for \eqn{\beta} (\eqn{k} columns), \eqn{\sigma^2} (1 column),
#'     \eqn{\theta} (1 column, if appropriate), \eqn{u} (\eqn{n} columns, not
#'     provided for log-normal model), \eqn{\log(t)} (\eqn{n} columns, simulated
#'     via data augmentation) and the logarithm of the adaptive variances (the
#'     number varies among models). The latter allows the user to evaluate if
#'     the adaptive variances have been stabilized.
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations (especially for the log-exponential power model).
#'
#' LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#'
#' @export
MCMC_LEP <- function(N,
                     thin,
                     burn,
                     Time,
                     Cens,
                     X,
                     beta0 = NULL,
                     sigma20 = NULL,
                     alpha0 = NULL,
                     prior = 2,
                     set = TRUE,
                     eps_l = 0.5,
                     eps_r = 0.5,
                     ar = 0.44) {

  # Sample starting values if not given
  if (is.null(beta0)) beta0 <- beta.sample(n = ncol(X))
  if (is.null(sigma20)) sigma20 <- sigma2.sample()
  if (is.null(alpha0)) alpha0 <- alpha.sample()

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

  beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
  beta[1, ] <- beta0
  sigma2 <- rep(0, times = N.aux + 1)
  sigma2[1] <- sigma20
  alpha <- rep(0, times = N.aux + 1)
  alpha[1] <- alpha0
  logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
  logt[1, ] <- log(Time)
  U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
  a <- ((abs(log(Time) - X %*% beta0)) / sqrt(sigma20)) ^ alpha0
  U0 <- -log(1 - stats::runif(n)) + a
  U[1, ] <- U0

  accept.beta <- rep(0, times = k)
  pbeta.aux <- rep(0, times = k)
  ls.beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
  accept.sigma2 <- 0
  psigma2.aux <- 0
  ls.sigma2 <- rep(0, times = N.aux + 1)
  accept.alpha <- 0
  palpha.aux <- 0
  ls.alpha <- rep(0, times = N.aux + 1)

  beta.aux <- beta[1, ]
  sigma2.aux <- sigma2[1]
  alpha.aux <- alpha[1]
  logt.aux <- logt[1, ]
  U.aux <- U[1, ]
  ls.beta.aux <- ls.beta[1, ]
  ls.sigma2.aux <- ls.sigma2[1]
  ls.alpha.aux <- ls.alpha[1]

  i_batch <- 0

  for (iter in 2:(N + 1)) {
    i_batch <- i_batch + 1

    for (ind.b in 1:k) {
      MH.beta <- MH_marginal_beta_j(omega2 = exp(ls.beta.aux[ind.b]),
                                    logt = logt.aux,
                                    X = X,
                                    sigma2 = sigma2.aux,
                                    alpha = alpha.aux,
                                    beta0 = beta.aux,
                                    j = ind.b)
      beta.aux[ind.b] <- MH.beta$beta[ind.b]
      if (MH.beta$ind == 1) {
        accept.beta[ind.b] <- accept.beta[ind.b] + 1
        pbeta.aux[ind.b] <- pbeta.aux[ind.b] + 1
      }
    }

    MH.sigma2 <- MH_marginal_sigma2(omega2 = exp(ls.sigma2.aux),
                                    logt = logt.aux,
                                    X = X,
                                    beta = beta.aux,
                                    alpha = alpha.aux,
                                    sigma20 = sigma2.aux,
                                    prior = prior)
    sigma2.aux <- MH.sigma2$sigma2
    if (MH.sigma2$ind == 1) {
      accept.sigma2 <- accept.sigma2 + 1
      psigma2.aux <- psigma2.aux + 1
    }

    MH.alpha <- MH_marginal_alpha(omega2 = exp(ls.alpha.aux),
                                  logt = logt.aux,
                                  X = X,
                                  beta = beta.aux,
                                  sigma2 = sigma2.aux,
                                  alpha0 = alpha.aux,
                                  prior = prior)
    alpha.aux <- MH.alpha$alpha
    if (MH.alpha$ind == 1) {
      accept.alpha <- accept.alpha + 1
      palpha.aux <- palpha.aux + 1
    }

    a <- ((abs(logt.aux - X %*% beta.aux)) / sqrt(sigma2.aux)) ^ alpha.aux
    U.aux <- -log(1 - stats::runif(n)) + a

    logt.aux <- logt.update.LEP(Time,
                                Cens,
                                X,
                                beta.aux,
                                sigma2.aux,
                                alpha.aux,
                                u = U.aux,
                                set,
                                eps_l,
                                eps_r)

    if (i_batch == 50) {
      pbeta.aux <- pbeta.aux / 50
      Pbeta.aux <- as.numeric(pbeta.aux < rep(ar, times = k))
      ls.beta.aux <- ls.beta.aux + ((-1)^Pbeta.aux) *
        min(0.01, 1 / sqrt(iter))
      psigma2.aux <- psigma2.aux / 50
      Psigma2.aux <- as.numeric(psigma2.aux < ar)
      ls.sigma2.aux <- ls.sigma2.aux + ((-1)^Psigma2.aux) *
        min(0.01, 1 / sqrt(iter))
      palpha.aux <- palpha.aux / 50
      Palpha.aux <- as.numeric(palpha.aux < ar)
      ls.alpha.aux <- ls.alpha.aux + ((-1)^Palpha.aux) *
        min(0.01, 1 / sqrt(iter))
      i_batch <- 0
      pbeta.aux <- rep(0, times = k)
      psigma2.aux <- 0
      palpha.aux <- 0
    }

    if (iter %% thin == 0) {
      beta[iter / thin + 1, ] <- beta.aux
      sigma2[iter / thin + 1] <- sigma2.aux
      alpha[iter / thin + 1] <- alpha.aux
      logt[iter / thin + 1, ] <- logt.aux
      U[iter / thin + 1, ] <- U.aux
      ls.beta[iter / thin + 1, ] <- ls.beta.aux
      ls.sigma2[iter / thin + 1] <- ls.sigma2.aux
      ls.alpha[iter / thin + 1] <- ls.alpha.aux
    }
    if ((iter - 1) %% 1e+05 == 0) {
      cat(paste("Iteration :", iter, "\n"))
    }
  }

  cat(paste("AR beta", 1:k, ":", round(accept.beta / N, 2), "\n"))
  cat(paste("AR sigma2 :", round(accept.sigma2 / N, 2), "\n"))
  cat(paste("AR alpha :", round(accept.alpha / N, 2), "\n"))

  chain <- cbind(beta, sigma2, alpha, U, logt, ls.beta, ls.sigma2, ls.alpha)

  beta.cols <- paste("beta.", seq(ncol(beta)), sep = "")
  alpha.cols <- paste("alpha.", seq(ncol(alpha)), sep = "")
  U.cols <- paste("U.", seq(ncol(U)), sep = "")
  logt.cols <- paste("logt.", seq(ncol(logt)), sep = "")
  ls.beta.cols <- paste("ls.beta.", seq(ncol(ls.beta)), sep = "")

  colnames(chain) <- c(beta.cols,
                       "sigma2",
                       alpha.cols,
                       U.cols,
                       logt.cols,
                       ls.beta.cols,
                       "ls.sigma2",
                       "ls.alpha")


  if (burn > 0) {
    burn.period <- 1:(burn / thin)
    chain <- chain [- burn.period, ]
  }

  return(chain)

}
#' @title Log-marginal likelihood estimator for the log-exponential power model
#' @inheritParams MCMC_LN
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=100 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations (especially for the log-exponential power model).
#'
#' LEP <- MCMC_LEP(N = 100, thin = 2, burn = 20, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' LEP.LML <- LML_LEP(thin = 2, Time = cancer[, 1], Cens = cancer[, 2],
#'                    X = cancer[, 3:11], chain = LEP)
#'
#' @export
LML_LEP <- function(thin,
                    Time,
                    Cens,
                    X,
                    chain,
                    prior = 2,
                    set = TRUE,
                    eps_l = 0.5,
                    eps_r = 0.5) {
  chain <- as.matrix(chain)
  n <- length(Time)
  N <- dim(chain)[1]
  k <- dim(X)[2]

  if (k > 1) {
    omega2.beta <- exp(apply(chain[, (2 * n + k + 3):(2 * n + 2 * k + 2)],
                             2, "median"))
  } else {
    omega2.beta <- exp(stats::median(chain[, 2 * n + 4]))
  }
  omega2.sigma2 <- exp(stats::median(chain[, 2 * n + 2 * k + 3]))
  omega2.alpha <- exp(stats::median(chain[, 2 * n + 2 * k + 4]))

  chain.nonadapt <- MCMC.LEP.NonAdapt(N = N * thin,
                                      thin = thin,
                                      Time,
                                      Cens,
                                      X,
                                      beta0 = as.vector(chain[N, 1:k]),
                                      sigma20 = chain[N, k + 1],
                                      alpha0 = chain[N, k + 2],
                                      prior,
                                      set,
                                      eps_l,
                                      eps_r,
                                      omega2.beta,
                                      omega2.sigma2,
                                      omega2.alpha)
  chain.nonadapt <- chain.nonadapt[-1, ]
  if (k > 1) {
    beta.star <- apply(chain.nonadapt[, 1:k], 2, "median")
  } else {
    beta.star <- stats::median(chain.nonadapt[, 1])
  }
  sigma2.star <- stats::median(chain.nonadapt[, k + 1])
  alpha.star <- stats::median(chain.nonadapt[, k + 2])

  # Log-LIKELIHOOD ORDINATE
  LL.ord <- log.lik.LEP(Time,
                        Cens,
                        X,
                        beta = beta.star,
                        sigma2 = sigma2.star,
                        alpha = alpha.star,
                        set,
                        eps_l,
                        eps_r)
  cat("Likelihood ordinate ready!\n")

  # PRIOR ORDINATE
  LP.ord <- prior_LEP(beta = beta.star,
                      sigma2 = sigma2.star,
                      alpha = alpha.star,
                      prior,
                      logs = TRUE)
  cat("Prior ordinate ready!\n")

  chain.sigma2 <- MCMCR.alpha.LEP(N = N * thin,
                                  thin = thin,
                                  Time, Cens,
                                  X,
                                  beta0 = as.vector(chain.nonadapt[N, 1:k]),
                                  sigma20 = chain.nonadapt[N, (k + 1)],
                                  alpha0 = alpha.star,
                                  logt0 = t(chain.nonadapt[N, (k + 3 + n) :
                                                             (2 * n + k + 2)]),
                                  u0 = t(chain.nonadapt[N,
                                                        (k + 3):(k + 2 + n)]),
                                  prior,
                                  set,
                                  eps_l,
                                  eps_r,
                                  omega2.beta,
                                  omega2.sigma2)
  cat("Reduced chain.sigma2 ready!\n")

  # POSTERIOR ORDINATE - alpha
  # USING AN ADAPTATION of the CHIB AND JELIAZKOV (2001) METHODOLOGY
  po1.alpha <- rep(0, times = N)
  po2.alpha <- rep(0, times = N)
  for (i in 1:N) {
    po1.alpha[i] <- alpha_alpha(alpha0 = as.numeric(chain.nonadapt[i,
                                                                   (k + 2)]),
                                alpha1 = alpha.star,
                                logt = as.vector(chain.nonadapt[i, (k + 3 + n) :
                                                                  (2 * n + k + 2)]),
                                X = X,
                                beta = as.vector(chain.nonadapt[i, 1:k]),
                                sigma2 = as.numeric(chain.nonadapt[i, (k + 1)]),
                                prior = prior) * stats::dnorm(x = alpha.star,
                                                              mean = as.numeric(chain.nonadapt[i, (k + 2)]),
                                                              sd = sqrt(omega2.alpha))
    alpha.aux <- stats::rnorm(n = 1, mean = alpha.star, sd = sqrt(omega2.alpha))
    po2.alpha[i] <- alpha_alpha(alpha0 = alpha.star, alpha1 = alpha.aux,
                                logt = as.vector(chain.sigma2[i, (k + 2 + n) :
                                                                (2 * n + k + 1)]),
                                X = X,
                                beta = as.vector(chain.sigma2[i + 1, 1:k]),
                                sigma2 = as.numeric(chain.sigma2[i + 1,
                                                                 (k + 1)]),
                                prior = prior)
  }
  PO.alpha <- mean(po1.alpha) / mean(po2.alpha)
  cat("Posterior ordinate alpha ready!\n")

  chain.beta <- MCMCR.sigma2.alpha.LEP(N = N * thin,
                                       thin = thin,
                                       Time,
                                       Cens,
                                       X,
                                       beta0 = as.vector(chain.nonadapt[N, 1:k]),
                                       sigma20 = sigma2.star,
                                       alpha0 = alpha.star,
                                       logt0 = t(chain.nonadapt[N, (k + 3 + n) :
                                                                  (2 * n + k + 2)]),
                                       u0 = t(chain.nonadapt[N, (k + 3) :
                                                               (k + 2 + n)]),
                                       prior,
                                       set,
                                       eps_l,
                                       eps_r,
                                       omega2.beta)
  cat("Reduced chain.beta ready\n!")

  # POSTERIOR ORDINATE - sigma2
  po1.sigma2 <- rep(0, times = N)
  po2.sigma2 <- rep(0, times = N)
  for (i in 1:N) {
    po1.sigma2[i] <- alpha_sigma2(sigma2_0 = as.numeric(chain.sigma2[i + 1,
                                                                     (k + 1)]),
                                  sigma2_1 = sigma2.star,
                                  logt = as.vector(chain.sigma2[i,
                                                                  (k + 2 + n) :
                                                                  (2 * n + k + 1)]),
                                  X = X,
                                  beta = as.vector(t(chain.sigma2[i + 1, 1:k])),
                                  alpha = alpha.star,
                                  prior = prior) *
      stats::dnorm(x = sigma2.star,
                   mean = as.numeric(chain.sigma2[i + 1, (k + 1)]),
                   sd = sqrt(omega2.sigma2))
    sigma2.aux <- stats::rnorm(n = 1,
                               mean = sigma2.star,
                               sd = sqrt(omega2.sigma2))
    po2.sigma2[i] <- alpha_sigma2(sigma2_0 = sigma2.star,
                                  sigma2_1 = sigma2.aux,
                                  logt = as.vector(chain.beta[i, (k + 1 + n):
                                                                (2 * n + k)]),
                                  X = X,
                                  beta = as.vector(chain.beta[i + 1, 1:k]),
                                  alpha = alpha.star,
                                  prior = prior)
  }
  PO.sigma2 <- mean(po1.sigma2) / mean(po2.sigma2)
  cat("Posterior ordinate sigma2 ready!\n")

  # POSTERIOR ORDINATE - beta
  chain.prev <- chain.beta
  PO.beta <- rep(0, times = k)

  for (j.beta in 0:(k - 1)) {
    print(j.beta)
    beta0 <- as.vector(chain.prev[N, 1:k])
    beta0[j.beta + 1] <- beta.star[j.beta + 1]
    chain.next <- MCMCR.betaJ.sigma2.alpha.LEP(N = N * thin,
                                               thin = thin,
                                               Time,
                                               Cens,
                                               X,
                                               beta0 = beta0,
                                               sigma20 = sigma2.star,
                                               alpha0 = alpha.star,
                                               logt0 = as.vector(chain.nonadapt[N,
                                                                                (k + 3 + n):
                                                                                  (2 * n + k + 2)]),
                                               u0 = as.vector(chain.nonadapt[N,
                                                                             (k + 3):
                                                                               (k + 2 + n)]),
                                               prior,
                                               set,
                                               eps_l,
                                               eps_r,
                                               omega2.beta,
                                               J = j.beta + 1)
    po1.beta <- rep(0, times = N)
    po2.beta <- rep(0, times = N)

    for (i in 1:N) {
      beta.0 <- as.vector(t(chain.prev[i + 1, 1:k]))
      beta.1 <- beta.0
      beta.1[j.beta + 1] <- beta.star[j.beta + 1]
      po1.beta[i] <- alpha_beta(beta_0 = beta.0,
                                beta_1 = beta.1,
                                logt = as.vector(chain.prev[i + 1,
                                                            (k + 1 + n) :
                                                              (2 * n + k)]),
                                X = X,
                                sigma2 = sigma2.star,
                                alpha = alpha.star) *
        stats::dnorm(x = beta.star[j.beta + 1],
                     mean = as.numeric(chain.prev[i + 1, j.beta + 1]),
                     sd = sqrt(omega2.beta[j.beta + 1]))
      betaj.aux <- stats::rnorm(n = 1, mean = beta.star[j.beta + 1],
                                sd = sqrt(omega2.beta[j.beta + 1]))
      beta.2 <- beta.star
      beta.2[j.beta + 1] <- betaj.aux
      po2.beta[i] <- alpha_beta(beta_0 = beta.star, beta_1 = beta.2,
                                logt = as.vector(chain.next[i + 1, (k + 1
                                                                    + n):
                                                              (2 * n + k)]),
                                X = X,
                                sigma2 = sigma2.star,
                                alpha = alpha.star)
    }
    PO.beta[j.beta + 1] <- mean(po1.beta) / mean(po2.beta)

    chain.prev <- chain.next
  }
  cat("Posterior ordinate beta ready!\n")

  # TAKING LOGARITHM
  LPO.alpha <- log(PO.alpha)
  LPO.sigma2 <- log(PO.sigma2)
  LPO.beta <- log(PO.beta)

  # MARGINAL LOG-LIKELIHOOD
  LML <- LL.ord + LP.ord - LPO.alpha - LPO.sigma2 - sum(LPO.beta)

  return(list(LL.ord = LL.ord,
              LP.ord = LP.ord,
              LPO.alpha = LPO.alpha,
              LPO.sigma2 = LPO.sigma2,
              LPO.beta = sum(LPO.beta),
              LML = LML))
}

#' @title Deviance information criterion for the log-exponential power model
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
#' # estimations (especially for the log-exponential power model).
#'
#' LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' LEP.DIC <- DIC_LEP(Time = cancer[, 1], Cens = cancer[, 2],
#'                    X = cancer[, 3:11], chain = LEP)
#'
#' @export
DIC_LEP <- function(Time,
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
    LL[iter] <- log.lik.LEP(Time,
                            Cens,
                            X,
                            beta = as.vector(chain[iter, 1:k]),
                            sigma2 = chain[iter, k + 1],
                            alpha = chain[iter, k + 2],
                            set,
                            eps_l,
                            eps_r)
  }

  aux <- apply(chain[, 1:(k + 2)], 2, "median")
  pd <- -2 * mean(LL) + 2 * log.lik.LEP(Time,
                                        Cens,
                                        X,
                                        beta = aux[1:k],
                                        sigma2 = aux[k + 1],
                                        alpha = aux[k + 2],
                                        set,
                                        eps_l,
                                        eps_r)
  pd.aux <- k + 2

  DIC <- -2 * mean(LL) + pd

  cat(paste("Effective number of parameters :", round(pd, 2), "\n"))
  cat(paste("Actual number of parameters :", pd.aux), "\n")
  return(DIC)
}

#' @title Case deletion analysis for the log-exponential power model
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
#' # estimations (especially for the log-exponential power model).
#'
#' LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' LEP.CD <- CaseDeletion_LEP(Time = cancer[, 1], Cens = cancer[, 2],
#'                            X = cancer[, 3:11], chain = LEP)
#'
#' @export
CaseDeletion_LEP <- function(Time,
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
      aux2[ITER] <- log.lik.LEP(Time[s],
                                Cens[s],
                                X[s, ],
                                beta = as.vector(chain[ITER, 1:k]),
                                sigma2 = chain[ITER, k + 1],
                                alpha = chain[ITER, k + 2],
                                set,
                                eps_l,
                                eps_r)
      aux1[ITER] <- exp(-aux2[ITER])
    }
    logCPO[s] <- -log(mean(aux1))
    KL.aux[s] <- mean(aux2)
    if (KL.aux[s] - logCPO[s] < 0) {
      cat(paste("Numerical problems for observation:", s, "\n"))
    }
  }
  KL <- abs(KL.aux - logCPO)
  CALIBRATION <- 0.5 * (1 + sqrt(1 - exp(-2 * KL)))
  return(cbind(logCPO, KL, CALIBRATION))
}


#' @title Outlier detection for observation for the log-exponential power model
#' @description This returns a unique number corresponding to the Bayes Factor
#'     associated to the test \eqn{M_0: \Lambda_{obs} = \lambda_{ref}} versus
#'     \eqn{M_1: \Lambda_{obs}\neq \lambda_{ref}} (with all other
#'     \eqn{\Lambda_j,\neq obs} free). The value of \eqn{\lambda_{ref}} is
#'     required as input. The user should expect long running times for the
#'     log-Student’s t model, in which case a reduced chain given
#'     \eqn{\Lambda_{obs} = \lambda_{ref}} needs to be generated
#' @inheritParams MCMC_LN
#' @param burn Burn-in period
#' @param ref Reference value \eqn{u_{ref}}. Vallejos & Steel recommends this
#'     value be set to \eqn{1.6 +1_\alpha} for the LEP model.
#' @param ar Optimal acceptance rate for the adaptive Metropolis-Hastings
#'     updates
#' @param obs Indicates the number of the observation under analysis
#' @param chain MCMC chains generated by a BASSLINE MCMC function
#' @examples
#' library(BASSLINE)
#'
#' # Please note: N=1000 is not enough to reach convergence.
#' # This is only an illustration. Run longer chains for more accurate
#' # estimations (especially for the log-exponential power model).
#'
#' LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
#'                 Cens = cancer[, 2], X = cancer[, 3:11])
#' alpha <- mean(LEP[, 11])
#' uref <- 1.6 + 1 / alpha
#' LEP.Outlier <- BF_u_obs_LEP(N = 100, thin = 20, burn =1 , ref = uref,
#'                             obs = 1, Time = cancer[, 1], Cens = cancer[, 2],
#'                             cancer[, 3:11], chain = LEP)
#'
#' @export
BF_u_obs_LEP <- function(N,
                         thin,
                         burn,
                         ref,
                         obs,
                         Time,
                         Cens,
                         X,
                         chain,
                         prior = 2,
                         set = TRUE,
                         eps_l = 0.5,
                         eps_r = 0.5,
                         ar = 0.44) {
  chain <- as.matrix(chain)
  aux1 <- Post.u.obs.LEP(obs, ref, X, chain)
  aux2 <- CFP.obs.LEP(N,
                      thin,
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

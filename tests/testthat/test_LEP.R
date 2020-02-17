################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LEP returns expected number of rows when burn = 5,thin = 10",{

  N <- 100
  thin <- 10
  burn <- 20
  LEP <- MCMC_LEP(N = N, thin = thin, burn = burn, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
  expect_equal(nrow(LEP), N / thin + 1 - burn/thin)
})

################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################



test_that("II_alpha same result in C++ as in R",{
  set.seed(1)
  II.alpha <- function(alpha) {
      aux <- (1 / alpha) * sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
      return(aux)
  }
  expect_equivalent(II_alpha(0.5), II.alpha(0.5))
})

test_that("I_alpha same result in C++ as in R",{
  set.seed(1)
  I.alpha <- function(alpha) {
   aux <- sqrt(1 / alpha ^ 3) * sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) +
                                        (1 + digamma(1 + 1 / alpha))^ 2  - 1)
    return(aux)
  }
  expect_equivalent(I_alpha(0.5), I.alpha(0.5))
})

test_that("J_alpha same result in C++ as in R",{
  set.seed(123)
  J.alpha <- function(alpha, k) {
      aux <- ((alpha * (alpha - 1) * gamma(1 - 1 / alpha) /
                 gamma(1 / alpha)) ^ (k / 2)) * (1 / alpha) *
                    sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
      return(aux)
  }
  expect_equivalent(J_alpha(2, 2), J.alpha(2, 2))
})

test_that("d_texp same result in C++ as in R",{

  dtexp <- function(x, rate, trunc) {
    if (x >= trunc) {
      aux <- rate * exp(-rate * (x - trunc))
    } else {
      aux <- 0
    }
    return(aux)
  }
  expect_equal(d_texp(x = 2.5,trunc = 1.5), dtexp(2.5, 1, 1.5))
  expect_equal(d_texp(x = 1.5,trunc = 2.5), 0)
})

test_that("MH_marginal_sigma2 same in C++ as in R",{
  MH.marginal.sigma2 <- function(N = 1, omega2, logt, X, beta, alpha,
                                 sigma20, prior) {
    k <- length(beta)
    n <- length(logt)
    if (prior == 1) {
      p <- 1 + k / 2
    }
    if (prior == 2) {
      p <- 1
    }
    if (prior == 3) {
      p <- 1
    }
    sigma2 <- rep(0, times = N + 1)
    sigma2[1] <- sigma20
    ind <- rep(0, times = N + 1)
    for (i.s in 1:N) {
      y.aux <- stats::rnorm(n = 1, mean = sigma2[i.s], sd = sqrt(omega2))
      if (y.aux <= 0) {
        sigma2[i.s + 1] <- sigma2[i.s]
        ind[i.s + 1] <- 2
      } else {
        l1 <- - ((n / 2) + p) * log(y.aux / sigma2[i.s])
        l2 <- - sum(abs(logt - X %*% beta) ^ alpha) *
          (y.aux ^ (-alpha / 2) - sigma2[i.s] ^ (-alpha / 2))
        aux <- l1 + l2
        u.aux <- stats::runif(1, 0, 1)

        # THE FOLLOWING THREE LINES AVOID CRUSHES
        if (is.na(aux)) {
          sigma2[i.s + 1] <- sigma2[i.s]
          ind[i.s + 1] <- 0
          cat("NA.sigma2\n")
        }
        if (is.na(aux) == FALSE && aux == -Inf) {
          sigma2[i.s + 1] <- sigma2[i.s]
          ind[i.s + 1] <- 0
          cat("-Inf.sigma2\n")
        }
        if (is.na(aux) == FALSE && aux == Inf) {
          sigma2[i.s + 1] <- y.aux
          ind[i.s + 1] <- 1
          cat("Inf.sigma2\n")
        }

        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
          sigma2[i.s + 1] <- y.aux
          ind[i.s + 1] <- 1
        }
        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
          sigma2[i.s + 1] <- sigma2[i.s]
          ind[i.s + 1] <- 0
        }
      }
    }
    sigma2 <- sigma2[-1]
    ind <- ind[-1]
    list(sigma2 = sigma2, ind = ind)
  }

  # Test params
  omega2 <- 0.5
  logt <- c (2,3)
  X <- matrix(seq(4), ncol = 2)
  beta <- c(1,2)
  alpha <- 1.2
  sigma20 <- 0.5
  prior <- 2


  set.seed(123)
  # N will always equal 1 for MH
  R.result <- MH.marginal.sigma2(N = 1, omega2, logt,X,beta,alpha,sigma20,prior)

  set.seed(123)
  Cpp.result <- MH_marginal_sigma2(omega2, logt,X,beta,alpha,sigma20, prior)

  expect_equal(R.result,Cpp.result)
})

test_that("alpha_beta same in C++ as in R",{

  alpha.beta <- function(beta.0, beta.1, logt, X, sigma2, alpha) {
    l1 <- -sum(abs((logt - X %*% beta.1) / sqrt(sigma2)) ^ alpha)
    l2 <- -sum(abs((logt - X %*% beta.0) / sqrt(sigma2)) ^ alpha)
    aux <- min(1, exp(l1 - l2))

    return(aux)
  }

  R.result <- alpha.beta(c(1,2), c(2,3), c(0.2,0.3),
                         matrix(c(1,2,3,4), ncol = 2), 0.2, 0.5 )

  Cpp.result <- alpha_beta(c(1,2), c(2,3), c(0.2,0.3),
                         matrix(c(1,2,3,4), ncol = 2), 0.2, 0.5 )

  expect_equal(R.result, Cpp.result)
  })

test_that("alpha_sigma2 same in C++ as in R",{
  alpha.sigma2 <- function(sigma2.0, sigma2.1, logt, X, beta, alpha, prior) {
    if (sigma2.1 <= 0) {
      aux <- 0
    } else {
      k <- dim(X)[2]
      n <- dim(X)[1]
      if (prior == 1) {
        p <- 1 + k / 2
      }
      if (prior == 2) {
        p <- 1
      }
      if (prior == 3) {
        p <- 1
      }
      l1 <- - ((n / 2) + p) * log(sigma2.1 / sigma2.0)
      l2 <- - sum(abs(logt - X %*% (beta))^alpha) *
        (sigma2.1 ^ (-alpha / 2) - sigma2.0 ^ (-alpha / 2))
      aux <- min(1, exp(l1 + l2))
    }

    return(aux)
  }

  R.result <- alpha.sigma2(0.3, 0.5, c(1, 2), matrix(seq(4), ncol = 2), c(1,2),
                           0.4, 2)

  Cpp.result <- alpha_sigma2(0.3, 0.5, c(1, 2), matrix(seq(4), ncol = 2),
                             c(1,2), 0.4, 2)

  expect_equal(R.result, Cpp.result)

})

test_that("MH_marginal_alpha same in C++ as in R",{
  MH.marginal.alpha <- function(N = 1, omega2, logt, X, beta, sigma2,
                                alpha0, prior) {
    k <- dim(X)[2]
    n <- dim(X)[1]
    alpha <- rep(0, times = N + 1)
    alpha[1] <- alpha0
    ind <- rep(0, times = N + 1)

    for (i.a in 1:N) {
      y.aux <- stats::rnorm(n = 1, mean = alpha[i.a], sd = sqrt(omega2))
      if (y.aux >= 2 || y.aux <= 1) {
        alpha[i.a + 1] <- alpha[i.a]
        ind[i.a + 1] <- 0
      } else {
        l1 <- n * log(y.aux / gamma(1 / y.aux)) -
          sum(abs((logt - X %*% beta) / sqrt(sigma2)) ^ y.aux)
        l0 <- n * log(alpha[i.a] / gamma(1 / alpha[i.a])) -
          sum(abs((logt - X %*% beta) / sqrt(sigma2)) ^ alpha[i.a])
        aux <- l1 - l0 + log(prior_alpha(y.aux, k, prior) /
                               prior_alpha(alpha[i.a], k, prior))
        u.aux <- stats::runif(1, 0, 1)

        # THE FOLLOWING THREE LINES AVOID CRUSHES
        if (is.na(aux)) {
          alpha[i.a + 1] <- alpha[i.a]
          ind[i.a + 1] <- 0
          cat("NA.alpha\n")
        }
        if (is.na(aux) == FALSE && aux == -Inf) {
          alpha[i.a + 1] <- alpha[i.a]
          ind[i.a + 1] <- 0
          cat("-Inf.alpha\n")
        }
        if (is.na(aux) == FALSE && aux == Inf) {
          alpha[i.a + 1] <- y.aux
          ind[i.a + 1] <- 1
          cat("Inf.alpha\n")
        }

        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
          alpha[i.a + 1] <- y.aux
          ind[i.a + 1] <- 1
        }
        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
          alpha[i.a + 1] <- alpha[i.a]
          ind[i.a + 1] <- 0
        }
      }
    }
    alpha <- alpha[-1]
    ind <- ind[-1]
    list(alpha = alpha, ind = ind)
  }

  R.result <- MH.marginal.alpha(N = 1, 0.2, c(1,2), matrix(seq(4), ncol=2), c(2,1),
                                0.3, 0.4, 2)
  Cpp.result <- MH_marginal_alpha(0.2, c(1,2), matrix(seq(4), ncol=2),
                                  c(2,1), 0.3, 0.4, 2)
  expect_equal(R.result, Cpp.result)

})

test_that("MH_marginal_beta_j same in C++ as in R",{
  MH.marginal.beta.j <- function(N = 1, omega2, logt, X, sigma2, alpha,
                                 beta0, j) {
    k <- length(beta0)
    beta <- matrix(rep(0, times = k * (N + 1)), ncol = k)
    beta[1, ] <- beta0
    ind <- rep(0, times = N + 1)
    for (i.b in 1:N) {
      y.aux <- beta[i.b, ]
      y.aux[j] <- stats::rnorm(n = 1, mean = beta[i.b, j], sd = sqrt(omega2))
      aux <- -sum(abs((logt - X %*% y.aux) / sqrt(sigma2)) ^ alpha) +
        sum(abs((logt - X %*% beta[i.b, ]) / sqrt(sigma2)) ^ alpha)
      u.aux <- stats::runif(1, 0, 1)
      if (is.na(aux)) {
        beta[i.b + 1, ] <- beta[i.b, ]
        ind[i.b + 1] <- 0
        cat("NA.beta\n")
      }
      if (is.na(aux) == FALSE && aux == -Inf) {
        beta[i.b + 1, ] <- beta[i.b, ]
        ind[i.b + 1] <- 0
        cat("-Inf.beta\n")
      }
      if (is.na(aux) == FALSE && aux == Inf) {
        beta[i.b + 1, ] <- y.aux
        ind[i.b + 1] <- 1
        cat("Inf.beta\n")
      }
      if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
        beta[i.b + 1, ] <- y.aux
        ind[i.b + 1] <- 1
      }
      if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
        beta[i.b + 1, ] <- beta[i.b, ]
        ind[i.b + 1] <- 0
      }
    }
    beta <- beta[-1, ]
    ind <- ind[-1]
    list(beta = beta, ind = ind)
  }

  set.seed(123)
  R.result <- MH.marginal.beta.j(N = 1, 0.2, c(1, 2), matrix(seq(4), ncol = 2),
                                 0.3, 0.4, c(1, 2), 1)

  set.seed(123)
  Cpp.result <- MH_marginal_beta_j(0.2, c(1, 2), matrix(seq(4), ncol = 2),
                                   0.3, 0.4, c(1, 2), 1)

  expect_equal(R.result, Cpp.result)
})

test_that("alpha_alpha same in C++ as in R",{
  alpha.alpha <- function(alpha0, alpha1, logt, X, beta, sigma2, prior) {
    k <- dim(X)[2]
    n <- dim(X)[1]
    if (alpha1 < 1 || alpha1 > 2) {
      aux <- 0
    } else {
      l1 <- n * log(alpha1 / gamma(1 / alpha1)) -
        ((1 / sqrt(sigma2)) ^ alpha1) * sum(abs(logt - X %*%
                                                   beta) ^ alpha1)
      l0 <- n * log(alpha0 / gamma(1 / alpha0)) -
        ((1 / sqrt(sigma2)) ^ alpha0) * sum((abs(logt - (X) %*%
                                                   (beta))) ^ alpha0)
      aux <- min(1, exp(l1 - l0) * prior_alpha(alpha1, k, prior) /
                   prior_alpha(alpha0, k, prior))
    }
    return(aux)
  }

  R.result <- alpha.alpha(1.1, 1.2, c(1,2), matrix(seq(4), ncol = 2), c(2,1),
                          3, 2)

  Cpp.result <- alpha_alpha(1.1, 1.2, c(1,2), matrix(seq(4), ncol = 2), c(2,1),
                            3, 2)
  expect_equal(R.result, Cpp.result)
})

test_that("d_normp same in C++ as in R",{
  #### DENSITY FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
  #### (BASED ON library normalp).
  ### POSSIBLE
  dnormp <- function(x, mu = 0, sigmap = 1, p = 2, log = FALSE) {
    if (!is.numeric(x) || !is.numeric(mu) || !is.numeric(sigmap) || !is.numeric(p))
      stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1)
      stop("p must be at least equal to one")
    if (min(sigmap) <= 0)
      stop("sigmap must be positive")
    cost <- 2 * p ^ (1 / p) * gamma(1 + 1 / p) * sigmap
    expon1 <- (abs(x - mu)) ^ p
    expon2 <- p * sigmap ^ p
    dsty <- (1 / cost) * exp(-expon1 / expon2)
    if (log == TRUE)
      dsty <- log(dsty)
    dsty
  }

  R.result <- dnormp(c(1,2), c(0,0) , c(1,1), c(2,2), log = T)
  Cpp.result <- d_normp(c(1,2), c(0,0) , c(1,1), c(2,2), logs = T)
  expect_equal(R.result, Cpp.result)
})

test_that("pnorm same in C++ as in R",{
  pnormp <- function(q, mu = 0, sigmap = 1, p = 2,
                     lower.tail = TRUE, log.pr = FALSE) {
    if (!is.numeric(q) || !is.numeric(mu) || !is.numeric(sigmap) || !is.numeric(p))
      stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1)
      stop("p must be at least equal to one")
    if (min(sigmap) <= 0)
      stop("sigmap must be positive")
    z <- (q - mu) / sigmap
    zz <- abs(z) ^ p
    zp <- stats::pgamma(zz, shape = 1 / p, scale = p)
    lzp <- stats::pgamma(zz, shape = 1 / p, scale = p, log = TRUE)
    zp <- ifelse(z < 0, 0.5 - exp(lzp - log(2)), 0.5 + exp(lzp - log(2)))
    if (log.pr == TRUE)
      zp <- ifelse(z < 0, log(0.5 - exp(lzp - log(2))),
                   log(0.5 + exp(lzp - log(2))))
    zp
  }

  ## Test for when log.pr is FALSE

  R.result <- pnormp(c(1, 2), c(2, 3), c(1, 1), c(2, 1))
  Cpp.result <- p_normp(c(1, 2), c(2, 3), c(1, 1), c(2, 1))

  expect_equal(R.result, Cpp.result)


  ## Test for when log.pr is TRUE

  R.result <- pnormp(c(1, 2), c(2, 3), c(1, 1), c(2, 1), log.pr = TRUE)
  Cpp.result <- p_normp(c(1, 2), c(2, 3), c(1, 1), c(2, 1), log_pr = TRUE)
  expect_equal(R.result, Cpp.result)


})

test_that("log_lik_LEP same in C++ as in R",{

  log.lik.LEP <- function(Time, Cens, X, beta, sigma2, alpha, set, eps_l,
                          eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    alpha <- rep(alpha, times = n)
    SP <- as.vector(sqrt(sigma2) * (1 / alpha) ^ (1 / alpha))
    if (set == 1) {
      aux1 <- (I(Time > eps_l) * log(p_normp(log(Time + eps_r),
                                             mu = MEAN,
                                             sigmap = SP,
                                             p = alpha) -
                                       p_normp(log(abs(Time - eps_l)),
                                               mu = MEAN, sigmap = SP,
                                               p = alpha)) +
                 (1 - I(Time > eps_l)) * log(p_normp(log(Time + eps_r),
                                                     mu = MEAN,
                                                     sigmap = SP,
                                                     p = alpha) - 0))
      aux2 <- log(1 - p_normp(log(Time), mu = MEAN, sigmap = SP, p = alpha))
      aux <- ifelse(Cens == 1, aux1, aux2)
    }
    if (set == 0) {
      aux <- Cens * (d_normp(log(Time), mu = MEAN, sigmap = SP,
                             p = alpha, logs = TRUE) - log(Time)) + (1 - Cens) *
        log(1 - p_normp(log(Time), mu = MEAN,
                        sigmap = SP, p = alpha))
    }
    return(sum(aux))
  }

  # Set = 1
  R.result <- log.lik.LEP(c(1, 2), c(1, 0), matrix(seq(4), ncol=2), c(1, 2),
                          0.2, 1.2, 1, 0.3, 0.4)
  Cpp.result <- log_lik_LEP(c(1, 2), c(1, 0), matrix(seq(4), ncol=2), c(1, 2),
                            0.2, 1.2, 1, 0.3,0.4)
  expect_equal(R.result, Cpp.result)
})

# test_that("Post_u_obs_LEP same in C++ as in R",{
#   LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
#                   Cens = cancer[,2], X = cancer[,3:11])
#   alpha <- mean(LEP[,11])
#   uref <- 1.6 + 1 / alpha
#
#
#   Post.u.obs.LEP <- function(obs, ref, X, chain) {
#     N <- dim(chain)[1]
#     n <- dim(X)[1]
#     k <- dim(X)[2]
#     aux1 <- rep(0, times = N)
#     aux2 <- rep(0, times = N)
#
#     for (iter in 1:N) {
#       trunc.aux <- (abs(chain[iter, (obs + k + 2 + n)] - X[obs, ] %*%
#                           as.vector(chain[iter, 1:k])) /
#                       sqrt(chain[iter, k + 1])) ^ (chain[iter, k + 2])
#
#       aux1[iter] <- d_texp(x = ref, trunc = trunc.aux)
#     }
#     aux <- mean(aux1)
#     return(aux)
#   }
#
#   R.result <- Post.u.obs.LEP(obs = 1, ref = uref, X = cancer[,3:11],
#                              chain = LEP )
#   Cpp.result <- Post_u_obs_LEP(obs = 1, ref = uref, X = cancer[,3:11],
#                                chain = LEP)
#
#  expect_equal(R.result, Cpp.result)
# })


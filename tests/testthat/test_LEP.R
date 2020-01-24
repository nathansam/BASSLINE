################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LEP returns expected number of rows when burn = 5,thin = 10",{

  N <- 100
  thin <- 10
  burn <- 5
  LEP <- MCMC_LEP(N = N, thin = thin, burn = burn, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
  expect_equal(nrow(LEP), N / thin + 1 - burn)
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
  N <- 5
  omega2 <- 0.5
  logt <- c (2,3)
  X <- matrix(seq(4), ncol = 2)
  beta <- c(1,2)
  alpha <- 1.2
  sigma20 <- 0.5
  prior <- 2


  set.seed(123)
  R.result <- MH.marginal.sigma2(N, omega2, logt,X,beta,alpha,sigma20,prior)

  set.seed(123)
  Cpp.result <- MH_marginal_sigma2(N, omega2, logt,X,beta,alpha,sigma20, prior)

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
  Cpp.result <- alpha_sigma2(0.3, 0.5, c(1, 2), matrix(seq(4), ncol = 2), c(1,2),
                             0.4, 2)

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

  R.result <- MH.marginal.alpha(5, 0.2, c(1,2), matrix(seq(4), ncol=2), c(2,1),
                                0.3, 0.4, 2)
  Cpp.result <- MH_marginal_alpha(5, 0.2, c(1,2), matrix(seq(4), ncol=2),
                                  c(2,1), 0.3, 0.4, 2)
  expect_equal(R.result, Cpp.result)

})

# test_that("MH_marginal_beta_j same in C++ as in R",{
#   MH.marginal.beta.j <- function(N = 1, omega2, logt, X, sigma2, alpha,
#                                  beta0, j) {
#     k <- length(beta0)
#     beta <- matrix(rep(0, times = k * (N + 1)), ncol = k)
#     beta[1, ] <- beta0
#     ind <- rep(0, times = N + 1)
#     for (i.b in 1:N) {
#       y.aux <- beta[i.b, ]
#       y.aux[j] <- stats::rnorm(n = 1, mean = beta[i.b, j], sd = sqrt(omega2))
#       aux <- -sum(abs((logt - X %*% y.aux) / sqrt(sigma2)) ^ alpha) +
#         sum(abs((logt - X %*% beta[i.b, ]) / sqrt(sigma2)) ^ alpha)
#       u.aux <- stats::runif(1, 0, 1)
#       if (is.na(aux)) {
#         beta[i.b + 1, ] <- beta[i.b, ]
#         ind[i.b + 1] <- 0
#         cat("NA.beta\n")
#       }
#       if (is.na(aux) == FALSE && aux == -Inf) {
#         beta[i.b + 1, ] <- beta[i.b, ]
#         ind[i.b + 1] <- 0
#         cat("-Inf.beta\n")
#       }
#       if (is.na(aux) == FALSE && aux == Inf) {
#         beta[i.b + 1, ] <- y.aux
#         ind[i.b + 1] <- 1
#         cat("Inf.beta\n")
#       }
#       if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
#         beta[i.b + 1, ] <- y.aux
#         ind[i.b + 1] <- 1
#       }
#       if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
#         beta[i.b + 1, ] <- beta[i.b, ]
#         ind[i.b + 1] <- 0
#       }
#     }
#     beta <- beta[-1, ]
#     ind <- ind[-1]
#     list(beta = beta, ind = ind)
#   }
#
#   set.seed(123)
#   R.result <- MH.marginal.beta.j(5, 0.2, c(1, 2), matrix(seq(4), ncol = 2),
#                                  0.3, 0.4, c(1, 2), 1)
#
#   set.seed(123)
#   Cpp.result <- MH_marginal_beta_j(5, 0.2, c(1, 2), matrix(seq(4), ncol = 2),
#                                    0.3, 0.4, c(1, 2), 1)
#
#   expect_equal(R.result, Cpp.result)
# })

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

  R.result <- alpha.alpha(1.1, 1.2, c(1,2), matrix(seq(4), ncol = 2), c(2,1), 3, 2)

  Cpp.result <- alpha_alpha(1.1, 1.2, c(1,2), matrix(seq(4), ncol = 2), c(2,1), 3, 2)
  expect_equal(R.result, Cpp.result)
})



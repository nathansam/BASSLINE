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

################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################

test_that("expected value for prior.LEP when log = FALSE",{
  prior.val <- prior.LEP(0.5, 1, 1, 2, F)
  expect_equivalent(round(prior.val, 4), -1.3209)
})

test_that("expected value for prior.LEP when log = TRUE",{
  prior.val <- prior.LEP(0.5, 1, 1, 2, T)
  expect_equivalent(round(prior.val, 4), 0.2783)
})


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



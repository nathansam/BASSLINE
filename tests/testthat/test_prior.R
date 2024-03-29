test_that("r_GIG same result in C++ as in R", {
  set.seed(123)
  n <- 1
  r <- 4.5
  y <- stats::rnorm(n)
  y <- y^2
  y <- 1 + ((y - sqrt(y * (4 * r + y))) / (2 * r))
  u <- stats::runif(n)
  aux <- rep(0, times = n)
  for (i in 1:n) {
    if (u[i] <= 1 / (1 + y[i])) {
      aux[i] <- r / y[i]
    } else {
      aux[i] <- r * y[i]
    }
  }
  set.seed(123)
  expect_equivalent(aux, r_GIG(4.5))
})

test_that("prior_alpha same result in C++ as in R", {
  prior.alpha <- function(alpha, k, prior) {
    if (prior == 1) {
      aux <- J_alpha(alpha, k)
    }
    if (prior == 2) {
      aux <- II_alpha(alpha)
    }
    if (prior == 3) {
      aux <- I_alpha(alpha)
    }
    return(aux)
  }

  for (prior in 1:3) {
    expect_equal(prior.alpha(c(1, 2), 1, prior), prior_alpha(c(1, 2), 1, prior))
  }
})

test_that("prior_LEP same result in C++ as in R", {
  prior.LEP <- function(beta, sigma2, alpha, prior, log) {
    if (log == FALSE) {
      aux <- prior_LN(beta, sigma2, prior, logs = FALSE) *
        prior_alpha(alpha, length(beta), prior)
    }
    if (log == TRUE) {
      aux <- prior_LN(beta, sigma2, prior, logs = TRUE) +
        log(prior_alpha(alpha, length(beta), prior))
    }
    return(aux)
  }

  for (prior in 1:3) {
    expect_equal(
      prior.LEP(c(1, 2), 1, c(1, 2), prior, F),
      prior_LEP(c(1, 2), 1, c(1, 2), prior, F)
    )

    expect_equal(
      prior.LEP(c(1, 2), 1, c(1, 2), prior, T),
      prior_LEP(c(1, 2), 1, c(1, 2), prior, T)
    )
  }
})

test_that("Unexpected arguments for prior funcs return 0", {
  expect_equal(prior_nu_single(1, prior = 1), 0)
  expect_equal(prior_alpha(c(1, 2), 1, 4), numeric(0))
  expect_equal(prior_alpha_single(2, 1, 4), 0)
})

test_that("J_alpha same in C++ as in R", {
  J.alpha <- function(alpha, k) {
    aux <- ((alpha * (alpha - 1) * gamma(1 - 1 / alpha) / gamma(1 / alpha))^
      (k / 2)) *
      (1 / alpha) *
      sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
    return(aux)
  }

  expect_equal(J.alpha(1.1, 1), J_alpha_single(1.1, 1))
  expect_equal(J.alpha(c(1.1, 1.2), 1), J_alpha(c(1.1, 1.2), 1))
  expect_equal(J_alpha(1.1, 1), J_alpha_single(1.1, 1))
})

test_that("I_alpha same in C++ as in R", {
  I.alpha <- function(alpha) {
    aux <- sqrt(1 / alpha^3) * sqrt((1 + 1 / alpha) *
      trigamma(1 + 1 / alpha) +
      (1 + digamma(1 + 1 / alpha))^2 - 1)
    return(aux)
  }

  expect_equal(I.alpha(1.1), I_alpha_single(1.1))
  expect_equal(I.alpha(c(1.1, 1.2)), I_alpha(c(1.1, 1.2)))
  expect_equal(I_alpha(1.1), I_alpha_single(1.1))
})


test_that("II_alpha same in C++ as in R", {
  II.alpha <- function(alpha) {
    aux <- (1 / alpha) * sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
    return(aux)
  }
  expect_equal(II.alpha(1.1), II_alpha_single(1.1))
  expect_equal(II.alpha(c(1.1, 1.2)), II_alpha(c(1.1, 1.2)))
  expect_equal(II_alpha(1.1), II_alpha_single(1.1))
})

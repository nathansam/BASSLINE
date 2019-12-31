test_that("Prior_LN returns expect value for prior = 1",{

  beta <- c(4,1,4)
  sigma2 <- 0.1

  p <- 1 + length(beta) / 2
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4,1,4), sigma2 = 0.1, prior = 1,
                             logs = T), aux )

})

test_that("Prior_LN returns expected value for prior = 2",{

  beta <- c(4,1,4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4,1,4), sigma2 = 0.1, prior = 2,
                             logs = T), aux )

})

test_that("Prior_LN returns expected value for prior = 3 & logs = FALSE",{

  beta <- c(4,1,4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * sigma2
  expect_equivalent(prior_LN(beta = c(4,1,4), sigma2 = 0.1, prior = 3,
                             logs = F), aux )

})



test_that("IIJ_nu same result in C++ as in R",{
  nu <- seq(1,10)
  aux <- sqrt(nu/ (nu + 3)) * sqrt(trigamma(nu / 2) -trigamma((nu + 1) / 2) -
                                     (2 * (nu + 3))/(nu * (nu + 1) ^ 2))
  expect_equivalent(IIJ_nu(nu), aux)
})

test_that("r_GIG same result in C++ as in R",{
  set.seed(123)
  n <- 50
  r <- 4.5
  y <- stats::rnorm(n)
  y <- y^2
  y <- 1 + ((y - sqrt(y * (4 * r + y))) / (2 * r))
  u <- stats::runif(n)
  aux <- rep(0, times = n)
  for (i in 1 : n) {
    if (u[i] <= 1/(1 + y[i])) {
      aux[i] <- r/y[i]
    } else {
      aux[i] <- r * y[i]
    }
  }
  set.seed(123)
  expect_equivalent(aux, r_GIG(4.5, 50))
})


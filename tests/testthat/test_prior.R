test_that("prior.nu returns expected value for nu = 1 & prior = 2",{
  prior.val <- prior.nu(1,2)
  expect_equal(round(prior.val,4), 0.1914)
})

test_that("prior.alpha returns expected value for alpha = 1, k = 1, prior = 1",{
  prior.val <- prior.alpha(1,1,1)
  expect_equal(round(prior.val,4), 1.3209)
})

test_that("prior.alpha returns expected value for alpha = 1, k = 1, prior = 2",{
  prior.val <- prior.alpha(1,1,2)
  expect_equal(round(prior.val,4), 1.3209)
})

test_that("prior.alpha returns expected value for alpha = 1, k = 1, prior = 3",{
  prior.val <- prior.alpha(1,1,3)
  expect_equal(round(prior.val,4), 1.9031)
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


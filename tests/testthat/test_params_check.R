test_that("Check for when N is negative", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = -1, thin = 10, burn = 1,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 1
    ),
    "N should be an integer greater than zero.\n"
  )
})

test_that("Check for when N is zero", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 0, thin = 10, burn = 1,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1, set = 1,
      eps_l = 0.5, eps_r = 0.5
    ),
    "N should be an integer greater than zero.\n"
  )
})

test_that("Check for when N is not an integer", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 0.2, thin = 10, burn = 1,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "N should be an integer greater than zero.\n"
  )
})


test_that("Check for when thin is negative", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = -1, burn = 1,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "thin should be a integer > 2.\n"
  )
})

test_that("Check for when burn is negative", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = -10,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "burn should be a non-negative integer.\n"
  )
})


test_that("Check for when N  is less than burn", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 6, thin = 2, burn = 10,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "N must be greater than burn.\n"
  )
})

test_that("Check for when burn is not a multiple of thin", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 11,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "burn must a multiple of thin.\n"
  )
})

test_that("Check for when N is not a multiple of thin", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 101, thin = 10, burn = 10,
      Time = c(1, 2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "N must a multiple of thin.\n"
  )
})

test_that("Check for when Time is negative", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = c(1, -5), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "All values in Time should be non-negative.\n"
  )
})


test_that("Check for when Time is not of correct length", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = c(1, 5, 6), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 1, eps_r = 1
    ),
    "Time is not the correct length.\n"
  )
})


test_that("Check for when Cens is not 0/1", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = c(1, 1), Cens = c(2, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1, set = 1,
      eps_l = 0.5, eps_r = 0.5
    ),
    "Cens should be either 0 or 1 for each observation\n"
  )
})

test_that("Check for when Cens is not of correct length", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = seq(2), Cens = c(0, 1, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "Cens is not the correct length.\n"
  )
})


test_that("Check for when X is not a matrix", {
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = seq(5),
      Cens = c(0, 1, 0, 0, 1),
      X = c(0, 1, 0, 0, 1),
      beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1,
      set = 1, eps_l = 0.5, eps_r = 0.5
    ),
    "X is not a matrix.\n"
  )
})

test_that("Check for when beta0 is not of correct length", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = seq(2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2, 0.3),
      sigma20 = 1, prior = 1, set = 1,
      eps_l = 0.5, eps_r = 0.5
    ),
    "beta0 is not the correct length.\n"
  )
})

test_that("Check for when prior is not an expected value", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = seq(2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 4, set = 1,
      eps_l = 0.5, eps_r = 0.5
    ),
    "prior should be 1, 2 or 3. See documentation\n"
  )
})

test_that("Check for when set is not an expected value", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_error(
    MCMC.param.check(
      N = 100, thin = 10, burn = 10,
      Time = seq(2), Cens = c(0, 1),
      X = mat, beta0 = c(0.2, 0.2),
      sigma20 = 1, prior = 1, set = 4,
      eps_l = 0.5, eps_r = 0.5
    ),
    "set should be 1 or 2. See documentation\n"
  )
})

test_that("Check for all values are correct", {
  mat <- matrix(seq(4), nrow = 2)
  testthat::expect_equal(MCMC.param.check(
    N = 100, thin = 10, burn = 10,
    Time = seq(2), Cens = c(0, 1),
    X = mat, beta0 = c(0.2, 0.2),
    sigma20 = 1, prior = 2, set = 1,
    eps_l = 1, eps_r = 0.5
  ), NULL)
})

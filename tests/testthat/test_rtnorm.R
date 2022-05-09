test_that("rtnorm same as msm version for finite limits", {
  set.seed(123)
  Cpp.result <- rtnorm(
    n = 1, lower = -1, upper = 1,
    mu = 0, sd = 1
  )
  Cpp.result <- as.vector(Cpp.result)

  set.seed(123)
  R.result <- msm::rtnorm(
    n = 1, lower = -1,
    upper = 1, mean = 0,
    sd = 1
  )
  expect_equal(Cpp.result, R.result)
})


test_that("rtnorm same as msm version for inf upper", {
  set.seed(123)
  Cpp.result <- rtnorm(
    n = 1, lower = -1, upper = Inf,
    mu = 0, sd = 1
  )
  Cpp.result <- as.vector(Cpp.result)

  set.seed(123)
  R.result <- msm::rtnorm(
    n = 1, lower = -1,
    upper = Inf, mean = 0,
    sd = 1
  )
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm same as msm version for inf lower", {
  set.seed(123)
  Cpp.result <- rtnorm(
    n = 1, lower = -Inf, upper = 1,
    mu = 0, sd = 1
  )
  Cpp.result <- as.vector(Cpp.result)

  set.seed(123)
  R.result <- msm::rtnorm(
    n = 1, lower = -Inf,
    upper = 1, mean = 0,
    sd = 1
  )
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm with inf limits same as rnorm", {
  set.seed(123)
  Cpp.result <- rtnorm(
    n = 1, lower = -Inf, upper = Inf,
    mu = 0, sd = 1
  )
  Cpp.result <- as.vector(Cpp.result)

  set.seed(123)
  R.result <- rnorm(1, 0, 1)
  expect_equal(Cpp.result, R.result)
})

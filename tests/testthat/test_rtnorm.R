test_that("rtnorm same as rtruncnorm when all arguments are vectors", {
  set.seed(123)
  Cpp.result <- rtnorm(n = 2, lower = c(-1, -1), upper = c(1, 1),
                       mu = c(0, 0), sd = c(1, 1))
  Cpp.result <- as.vector(Cpp.result)
  
  set.seed(123)
  R.result <- truncnorm::rtruncnorm(n = 2, a = c(-1, -1),
                                    b = c(1, 1), mean = c(0, 0),
                                    sd = c(1, 1))
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm same as rtruncnorm when none of the arguments are vectors", {
  set.seed(123)
  Cpp.result <- rtnorm(n = 5, lower = -1, upper = Inf,
                       mu = 0, sd = 1)
  Cpp.result <- as.vector(Cpp.result)
  
  set.seed(123)
  R.result <- truncnorm::rtruncnorm(n = 5, a = -1,
                                    b = Inf, mean = 0,
                                    sd = 1)
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm same as rtruncnorm when SD isn't a vector", {

  #  Test when SD isn't a vector
  set.seed(123)
  Cpp.result <- rtnorm(n = 2, lower = c(-1, -1), upper = c(Inf, Inf),
                       mu = c(0, 0), sd = 1)
  Cpp.result <- as.vector(Cpp.result)
  
  set.seed(123)
  R.result <- truncnorm::rtruncnorm(n = 2, a = c(-1, -1),
                                    b = c(Inf, Inf), mean = c(0, 0), sd = 1)
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm same as rtruncnorm when lower/a isn't a vector", {
  # Test when lower / a isn't a vector
  set.seed(123)
  Cpp.result <- rtnorm(n = 2, lower = -1, upper = c(1, 1), mu = c(0, 0),
                       sd = c(1, 1))
  Cpp.result <- as.vector(Cpp.result)
  
  set.seed(123)
  R.result <- truncnorm::rtruncnorm(n = 2, a = -1,
                                    b = c(1, 1), mean = c(0, 0),
                                    sd = c(1, 1))
  expect_equal(Cpp.result, R.result)
})
  
test_that("rtnorm same as rtruncnorm when upper/b isn't a vector", {
  # Test when upper / b isn't a vector
  set.seed(123)
  Cpp.result <- rtnorm(n = 2, lower = c(-1, -1), upper = Inf, mu = c(0, 0),
                       sd = c(1, 1))
  Cpp.result <- as.vector(Cpp.result)
  set.seed(123)
  R.result <- truncnorm::rtruncnorm(n = 2, a = c(-1, -1),
                                    b = Inf, mean = c(0, 0),
                                    sd = c(1, 1))
  
  expect_equal(Cpp.result, R.result)
})

test_that("rtnorm returns same result when called from R as from C++", {
  R.result <- rtnorm(n = 3, lower = -Inf, upper = 1, mu = 0, sd = 1)
})
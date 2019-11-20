test_that("Prior_LN returns expect value for prior = 1",{

  beta <- c(4,1,4)
  sigma2 <- 0.1

  p <- 1 + length(beta) / 2
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4,1,4), sigma2 = 0.1, prior = 1,
                             logs = T), aux )

})

test_that("Prior_LN returns expect value for prior = 2",{

  beta <- c(4,1,4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4,1,4), sigma2 = 0.1, prior = 2,
                             logs = T), aux )

})

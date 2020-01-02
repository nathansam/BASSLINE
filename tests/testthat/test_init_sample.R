test_that("Check beta.sample returns expected values",{
  set.seed(123)
  betas <- round(beta.sample(4),2)
  expect_equal(betas, c(-0.56, -0.23, 1.56, 0.07))
})

test_that("Check sigma2.sample returns expected value",{
  set.seed(123)
  sigma2 <- round(sigma2.sample(),2)
  expect_equal(sigma2, 0.45)
})

test_that("Check nu.sample returns expected value",{
  set.seed(123)
  nu <- round(nu.sample(),2)
  expect_equal(nu, 0.45)
})

test_that("Check alpha.sample returns expected value",{
  set.seed(123)
  alpha <- round(alpha.sample(),2)
  expect_equal(alpha, 1.29)
})

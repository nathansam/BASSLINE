test_that("Expected value for prior_LST when logs = FALSE",{
  expect.val <- round( prior.LST(1,1,1,2,F),4)
  expect_equal(expect.val,-0.1914)
})

test_that("Expected value for prior_LST when logs = TRUE",{
  expect.val <- round( prior.LST(1,1,1,2,T),4)
  expect_equal(expect.val, -1.6536)
})

test_that("Expected value for MH.nu.LST",{
  set.seed(123)
  val <- MH.nu.LST(N = 1, omega2 = 1, beta = 1, lambda = 1, nu0 = 1,
                   prior = 2) $nu

  expect_equal(round(val, 4), 0.4395)
})

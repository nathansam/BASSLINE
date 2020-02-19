################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LLAP returns expected num of rows when burn = 5 ,thin = 10", {

  N <- 100
  thin <- 10
  burn <- 20
  LLAP <- MCMC_LLAP(N = N, thin = thin, burn = burn,  Time = cancer[, 1],
                    Cens = cancer[, 2], X = cancer[, 3:11])
  expect_equal(nrow(LLAP), N / thin + 1 - burn / thin)
})

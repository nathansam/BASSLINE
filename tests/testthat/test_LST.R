#test_that("MCMC_LST Converges",{
#  set.seed(123)
#  LST.chain <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
#                      Cens = cancer[,2], X = cancer[,3:4])
#
#  testthat::expect_true(abs(LN.chain[10,5] - LN.chain[11,5]) < 0.01)
#})


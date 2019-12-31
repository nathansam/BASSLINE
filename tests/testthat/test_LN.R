test_that("MCMC_LN Converges",{
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")

  LN.chain <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = c(5,10),
                      Cens = c(1,0), X = matrix(c(1,2,3,4), nrow = 2))

  testthat::expect_true(abs(LN.chain[10,5] - LN.chain[11,5]) < 0.01)
})


test_that("LML_LN Returns Expected Result",{
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")
  LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                Cens = cancer[,2], X = cancer[,3:11])
  LN.LML <- LML_LN(thin = 20, Time = cancer[,1], Cens = cancer[,2],
                   X = cancer[,3:11], chain = LN)
  testthat::expect_equal(as.numeric(LN.LML), c( -716.97105391, -0.06768493,
                                               0.77551463, 14.37110198,
                                               -732.18535545))
})

test_that("DIC_LN Returns Expected Result",{
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")
  LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                Cens = cancer[,2], X = cancer[,3:11])
  LN.DIC <- DIC_LN(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                   chain = LN)
  testthat::expect_equal(round(LN.DIC, 4), 1445.1216)
})

test_that("CaseDeletion_LN Returns Expected Result",{
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")
  LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                Cens = cancer[,2], X = cancer[,3:11])
  LN.CD <- CaseDeletion_LN(Time = cancer[,1], Cens = cancer[,2],
                           X = cancer[,3:11], chain = LN)
  means <- round(c(mean(LN.CD[,1]), mean(LN.CD[,2]), mean(LN.CD[,3])), 4)
  testthat::expect_equal(means, c(-5.2956, 0.0418, 0.6029))
})

test_that("MCMC_LN Converges",{
  set.seed(123)
  LN.chain <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = c(5,10),
                      Cens = c(1,0), X = matrix(c(1,2,3,4), nrow = 2))

  testthat::expect_true(abs(LN.chain[10,5] - LN.chain[11,5]) < 0.01)
})


test_that("LML_LN Returns Expected Result",{
  if(.Machine$sizeof.pointer == 8){
    set.seed(123)
    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
    LN.LML <- LML_LN(thin = 20, Time = cancer[,1], Cens = cancer[,2],
                     X = cancer[,3:11], chain = LN)
    testthat::expect_equal(round(as.numeric(LN.LML),2), c( -717.22, -0.27,
                                                           0.75, 13.94,
                                                           -732.17))
    }
})

test_that("DIC_LN Returns Expected Result",{
  if(.Machine$sizeof.pointer == 8){
    set.seed(123)
    
    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
    LN.DIC <- DIC_LN(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                     chain = LN)
    testthat::expect_equal(round(LN.DIC, 4), 1449.8554)
    }
})

test_that("CaseDeletion_LN Returns Expected Result",{
  if(.Machine$sizeof.pointer == 8){
    set.seed(123)
    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
    LN.CD <- CaseDeletion_LN(Time = cancer[,1], Cens = cancer[,2],
                             X = cancer[,3:11], chain = LN)
    means <- round(c(mean(LN.CD[,1]), mean(LN.CD[,2]), mean(LN.CD[,3])), 4)
    testthat::expect_equal(means, c(-5.312, 0.0487, 0.6068))}
})

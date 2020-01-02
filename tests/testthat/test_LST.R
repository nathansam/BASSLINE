################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LST returns expected number of rows when burn = 0",{
  N <- 1000
  thin <- 20
  LST <- MCMC_LST(N = N, thin = thin, burn = 0 , Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])
  expect_equal(nrow(LST), N / thin + 1)
})

################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################

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

test_that("Expected value for alpha.nu when nu1 > 0",{
  if(.Machine$sizeof.pointer == 8){

    set.seed(123)
    val <- alpha.nu(1, nu1 =2, 1, 1, 2)

    expect_equal(round(val, 4), 0.6657)
  }
})

test_that("Expected value for alpha.nu when nu1 == 0",{
  set.seed(123)
  val <- alpha.nu(1, nu1 =0, 1, 1, 2)
  expect_equal(round(val, 4), 0)
})

test_that("log.lik.LST returns expected value for set = 0",{
  if(.Machine$sizeof.pointer == 8){
    set.seed(123)
    lik <- log.lik.LST(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                       beta = rep(0,9), sigma2 = 1, nu = 1, set = 0 ,
                       eps_l = 0.5, eps_r = 0.5)

    expect_equal(round(lik,4), -1043.1982)
  }
})

test_that("log.lik.LST returns expected values for set = 1",{
  if(.Machine$sizeof.pointer == 8){
    set.seed(123)
    lik <- log.lik.LST(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                       beta = rep(0,9), sigma2 = 1, nu = 1, set = 1 ,
                       eps_l = 0.5, eps_r = 0.5)

    expect_equal(round(lik,4), -1043.0675)
  }
})

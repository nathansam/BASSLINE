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

test_that("prior_LST is the same in C++ as in R",{
  prior.LST.R <- function(beta, sigma2, nu, prior, log) {
    if (log == FALSE) {
      aux <- prior_LN(beta, sigma2, prior, logs = FALSE) *
        prior_nu(nu, prior)
    }
    if (log == TRUE) {
      aux <- prior_LN(beta, sigma2, prior, logs = TRUE) +
        log(prior_nu(nu, prior))
    }
    return(aux)
  }

  expect_equal(prior.LST.R(c(1,2),1,1,2,T), prior_LST(c(1,2),1,1,2,T))
  expect_equal(prior.LST.R(c(1,2),1,1,2,F), prior_LST(c(1,2),1,1,2,F))
})



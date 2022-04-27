################################################################################
################################ USER FUNCTIONS ################################
################################################################################


test_that("MCMC_LN returns expected number of rows when burn = 5,thin = 10", {

  N <- 100
  thin <- 10
  burn <- 10
  LN <- MCMC_LN(N = N, thin = thin, burn = burn, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
  expect_equal(nrow(LN), N / thin + 1 - (burn / thin))
})




################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################

test_that("Prior_LN returns expect values", {


  prior.LN<-function(beta,sigma2,prior,log){
    k <- length(beta)
    if(prior == 1) p <- 1 + k / 2
    if(prior == 2) p <- 1
    if(prior==3) p <- 1
    if(log == FALSE) aux <- sigma2 ^ ( - p)
    if(log == TRUE) aux <- -p * log(sigma2)
    return(aux)
  }

  betas <- c(4, 1, 4)
  sigma2 <- 0.1

  for (prior in 1:3){
    for (logs in c(T, F)){
      expect_equivalent(prior_LN(betas, sigma2, prior, logs),
                        prior.LN(betas, sigma2, prior, logs ))
    }
  }
})

test_that("Prior_LN returns expected value for prior = 2", {

  beta <- c(4, 1, 4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4, 1, 4), sigma2 = 0.1, prior = 2,
                             logs = T), aux)

})


test_that("log.lik.LN returns expected value for set = 0", {
  if (.Machine$sizeof.pointer == 8) {
    set.seed(123)
    lik <- log.lik.LN(Time = cancer[, 1], Cens = cancer[, 2],
                      X = cancer[, 3:11], beta = rep(0, 9), sigma2 = 1,
                      set = 0, eps_l = 0.5, eps_r = 0.5)

    expect_equal(round(lik, 4), -1927.0135)
  }
})

test_that("log.lik.LN returns expected values for set = 1", {
  if (.Machine$sizeof.pointer == 8) {
    set.seed(123)
    lik <- log.lik.LN(Time = cancer[, 1], Cens = cancer[, 2],
                      X = cancer[, 3:11], beta = rep(0, 9), sigma2 = 1,
                      set = 1, eps_l = 0.5, eps_r = 0.5)

    expect_equal(round(lik, 4), -1926.7085)
  }
})



test_that("logt_update_SMLN same in C++ as in R",{
  Time <- cancer[, 1]
  Cens <- cancer[, 2]
  X <- cancer[, 4:11]
  beta <- beta.sample(8)
  sigma2 <- 1.5; set <- TRUE; eps_l <- 0.5; eps_r <- 0.5


  set.seed(123)
  result.R <- logt.update.SMLN(Time, Cens, X, beta , sigma2, set, eps_l,
                               eps_r)
  set.seed(123)
  result.Cpp <- logt_update_SMLN(Time, Cens, X, beta, sigma2, set, eps_l,
                                 eps_r)
  
  diff.within.tolerance <- mean(result.R - result.Cpp) < 0.03

  testthat::expect_equal(diff.within.tolerance, TRUE)

})



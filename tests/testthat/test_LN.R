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

test_that("LML_LN Returns Expected Result", {
  if (.Machine$sizeof.pointer == 8) {

    set.seed(123)
    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
    LN.LML <- LML_LN(thin = 20, Time = cancer[, 1], Cens = cancer[, 2],
                     X = cancer[, 3:11], chain = LN)
    testthat::expect_equal(round(as.numeric(LN.LML), 2), c(-715.29, -0.20,
                                                           0.85, 15.54,
                                                           -731.89))
    }
})

test_that("DIC_LN Returns Expected Result", {
  if (.Machine$sizeof.pointer == 8) {
    set.seed(123)

    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
    LN.DIC <- DIC_LN(Time = cancer[, 1], Cens = cancer[, 2], X = cancer[, 3:11],
                     chain = LN)
    testthat::expect_equal(round(LN.DIC, 4), 1450.0062)
    }
})

test_that("CaseDeletion_LN Returns Expected Result", {

  if (.Machine$sizeof.pointer == 8) {
    set.seed(123)
    LN <- MCMC_LN(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
    LN.CD <- CaseDeletion_LN(Time = cancer[, 1], Cens = cancer[, 2],
                             X = cancer[, 3:11], chain = LN)
    means <- round(c(mean(LN.CD[, 1]), mean(LN.CD[, 2]), mean(LN.CD[, 3])), 4)
    testthat::expect_equal(means, c(-5.2990, 0.0425, 0.6047))
  }
})

################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################

test_that("Prior_LN returns expect value for prior = 1", {

  beta <- c(4, 1, 4)
  sigma2 <- 0.1

  p <- 1 + length(beta) / 2
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4, 1, 4), sigma2 = 0.1, prior = 1,
                             logs = T), aux)

})

test_that("Prior_LN returns expected value for prior = 2", {

  beta <- c(4, 1, 4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * log(sigma2)
  expect_equivalent(prior_LN(beta = c(4, 1, 4), sigma2 = 0.1, prior = 2,
                             logs = T), aux)

})

test_that("Prior_LN returns expected value for prior = 3 & logs = FALSE", {

  beta <- c(4, 1, 4)
  sigma2 <- 0.1
  p <- 1
  aux <- -p * sigma2
  expect_equivalent(prior_LN(beta = c(4, 1, 4), sigma2 = 0.1, prior = 3,
                             logs = F), aux)
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

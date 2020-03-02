################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LST returns expected number of rows when burn = 0", {
  N <- 1000
  thin <- 20
  LST <- MCMC_LST(N = N, thin = thin, burn = 0, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
  expect_equal(nrow(LST), N / thin + 1)
})

test_that("LML_LST Returns Expected Result", {
  if (.Machine$sizeof.pointer == 8) {

    set.seed(123)
    LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
    LST.LML <- LML_LST(thin = 20, Time = cancer[, 1], Cens = cancer[, 2],
                     X = cancer[, 3:11], chain = LST)
    #testthat::expect_equal(round(as.numeric(LST.LML),2), c( -714.06, -2.16,
    #                                                       -1.62, 1.15,
    #                                                       14.52, -730.26))
  }
})


test_that("DIC_LST Returns Expected Result", {
  if (.Machine$sizeof.pointer == 8) {

    set.seed(123)
    LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                    Cens = cancer[, 2], X = cancer[, 3:11])
    LST.DIC <- DIC_LST(Time = cancer[, 1], Cens = cancer[, 2],
                       X = cancer[, 3:11], chain = LST)
    #testthat::expect_equal(round(LST.DIC,2), 1445.54)
  }
})


################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################

test_that("Expected value for alpha_nu when nu1 > 0", {
  if (.Machine$sizeof.pointer == 8) {

    set.seed(123)
    val <- alpha_nu(1, nu1 = 2, 1, 2)

    expect_equal(round(val, 4), 0.6657)
  }
})

test_that("Expected value for alpha_nu when nu1 == 0", {
  set.seed(123)
  val <- alpha_nu(1, nu1 = 0, 1, 2)
  expect_equal(val, 0)
})

test_that("log_lik_LST same in C++ as in R", {

  log.lik.LST <- function(Time, Cens, X, beta, sigma2, nu, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    nu <- rep(nu, times = n)
    if (set == 1) {
      aux <- Cens * (I(Time > eps_l) *
                       log(stats::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2),
                                     df = nu) -
                             stats::pt((log(abs(Time - eps_l)) - MEAN) / sqrt(sigma2),
                                       df = nu)) +
                       (1 - I(Time > eps_l)) *
                       log(stats::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2),
                                     df = nu) - 0)) +
        (1 - Cens) *
        log(1 - stats::pt((log(Time) - MEAN) / sqrt(sigma2),
                          df = nu))
    }
    if (set == 0) {
      aux <- Cens * (stats::dt((log(Time) - MEAN) / sqrt(sigma2),
                               df = nu, log = TRUE) -
                       log(sqrt(sigma2) * Time)) +
        (1 - Cens) * log(1 - stats::pt((log(Time) - MEAN) / sqrt(sigma2),
                                       df = nu))
    }
    return(sum(aux))
  }

  if (.Machine$sizeof.pointer == 8) {

    for (set in 0:1) {

      set.seed(123)
      R.result <- log.lik.LST(Time = cancer[, 1], Cens = cancer[, 2],
                              X = cancer[, 3:11], beta = seq(9),
                              sigma2 = 1, nu = 1, set = set,
                              eps_l = 0.5, eps_r = 0.5)
      set.seed(123)
      Cpp.result <- log_lik_LST(Time = cancer[, 1], Cens = cancer[, 2],
                                X = cancer[, 3:11], beta = seq(9),
                                sigma2 = 1, nu = 1, set = set,
                                eps_l = 0.5, eps_r = 0.5)

      expect_equal(R.result, Cpp.result)
    }
  }
})

test_that("prior_LST is the same in C++ as in R", {
  prior.LST.R <- function(beta, sigma2, nu, prior, log) {
    if (log == FALSE) {
      aux <- prior_LN(beta, sigma2, prior, logs = FALSE) *
        prior_nu_single(nu, prior)
    }
    if (log == TRUE) {
      aux <- prior_LN(beta, sigma2, prior, logs = TRUE) +
        log(prior_nu_single(nu, prior))
    }
    return(aux)
  }


  expect_equal(prior.LST.R(c(1, 2), 1, 1, 2, T),
               prior_LST(c(1, 2), 1, 1, 2, T))
  expect_equal(prior.LST.R(c(1, 2), 1, 1, 2, F),
               prior_LST(c(1, 2), 1, 1, 2, F))
})

test_that("alpha_nu same in C++ as in R", {
  alpha.nu <- function(nu0, nu1, lambda, prior) {
    if (nu1 <= 0) {
      aux <- 0
    }
    if (nu1 > 0) {
      aux <- min(1, exp(sum(stats::dgamma(lambda, shape = nu1 / 2,
                                          rate = nu1 / 2, log = TRUE))
                        - sum(stats::dgamma(lambda,
                                            shape = nu0 / 2,
                                            rate = nu0 / 2,
                                            log = TRUE))) *
                   (prior_nu_single(nu1, prior) / prior_nu_single(nu0, prior)))
    }
    return(aux)
  }

  for (prior in 1:3) {
    expect_equal(alpha.nu(1, 2, c(1, 2), prior), alpha_nu(1, 2, c(1, 2), prior))

  }



})

test_that("MH_nu_LST same in C++ as in R", {
  ### METROPOLIS-HASTINMCMC UPDATE OF NU
  ### (REQUIRED FOR SEVERAL .LST FUNCTIONS)
  MH.nu.LST <- function(N = 1, omega2, beta, lambda, nu0, prior) {
    k <- length(beta)
    ind <- rep(0, N + 1)
    nu <- rep(0, times = N + 1)
    nu[1] <- nu0

    for (j.nu in 1:N) {
      y <- stats::rnorm(1, mean = nu[j.nu], sd = sqrt(omega2))
      ind.aux <- I(y >= 0)
      y <- abs(y)
      u.aux <- stats::runif(1, min = 0, max = 1)

      log.aux <- sum(stats::dgamma(lambda, shape = y / 2, rate = y / 2,
                                   log = TRUE)) -
        sum(stats::dgamma(lambda,
                          shape = nu[j.nu] / 2,
                          rate = nu[j.nu] / 2,
                          log = TRUE)) +
        log(prior_nu_single(y, prior) /
              prior_nu_single(nu[j.nu], prior))

      aux <- I(log(u.aux) < log.aux)
      aux <- as.numeric(aux) * as.numeric(ind.aux)
      nu[j.nu + 1] <- aux * y + (1 - aux) * nu[j.nu]
      ind[j.nu + 1] <- as.numeric(aux)
    }
    nu <- nu[-1]
    ind <- ind[-1]
    list(nu = nu, ind = ind)
  }

  for (prior in 2) {

    set.seed(123)
    MH.nu.LST.R <- MH.nu.LST(N = 1, omega2 = 0.5, beta = c(1, 2),
                             lambda = c(3, 4), nu0 = 0.2, prior = prior)

    set.seed(123)
    MH.nu.LST.Cpp <- MH_nu_LST(omega2 = 0.5, beta = c(1, 2),
                               lambda = c(3, 4), nu0 = 0.2, prior = prior)

    expect_equal(MH.nu.LST.R, MH.nu.LST.Cpp)
  }

})

test_that("Post_lambda_obs_LST same in C++ as in R", {

  Post.lambda.obs.LST <- function(obs, ref, X, chain) {
    N <- dim(chain)[1]
    n <- dim(X)[1]
    k <- dim(X)[2]
    aux1 <- rep(0, times = N)
    aux2 <- rep(0, times = N)

    for (iter in 1:N) {
      aux1[iter] <- (((chain[iter, (obs + k + 2 + n)] - X[obs, ] %*%
                         as.vector(chain[iter, (1:k)]))^2) / (chain[iter,
                                                                    (k + 1)])) +
        chain[iter, (k + 2)]
      aux2[iter] <- stats::dgamma(x = ref,
                                  shape = (chain[iter, (k + 2)] + 1) / 2,
                                  rate = aux1[iter] / 2)
    }
    aux <- mean(aux2)
    return(aux)
  }

  LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])


  R.result <- Post.lambda.obs.LST(1, 1, cancer[, 3:11], LST)
  Cpp.result <- Post_lambda_obs_LST(1, 1, cancer[, 3:11], LST)
  expect_equal(R.result, Cpp.result)
})

test_that("Prior_nu_single same in C++ as in R", {

  IIJ.nu <- function(nu, prior)
  {
    if (prior == 2){
    aux <- sqrt(nu / (nu + 3)) *
      sqrt(trigamma(nu / 2) - trigamma((nu + 1) / 2) - (2 * (nu + 3)) /
        (nu * (nu + 1) ^ 2))
    return(aux)
    }
  }

  nu <- runif(10, 0, 10)

  prior <- 2

  for (i in 1:10){
    aux <- IIJ.nu(nu[i], prior)
    expect_equal(aux, prior_nu_single(nu[i], prior))
  }
})

test_that("prior_nu_single returns expected value for nu = 1 & prior = 2", {
  prior.val <- prior_nu_single(1, 2)
  expect_equal(round(prior.val, 4), 0.5679)
})



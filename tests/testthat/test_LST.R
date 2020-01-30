################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LST returns expected number of rows when burn = 0",{
  N <- 1000
  thin <- 20
  LST <- MCMC_LST(N = N, thin = thin, burn = 0 , Time = cancer[, 1],
                  Cens = cancer[, 2], X = cancer[, 3:11])
  expect_equal(nrow(LST), N / thin + 1)
})

################################################################################
############################## INTERNAL FUNCTIONS ##############################
################################################################################



test_that("Expected value for MH_nu_LST",{
  set.seed(123)
  val <- MH_nu_LST(omega2 = 1, beta = 1, lambda = 1, nu0 = 1,
                   prior = 2) $nu

  expect_equal(round(val, 4), 0.4395)
})

test_that("Expected value for alpha_nu when nu1 > 0",{
  if(.Machine$sizeof.pointer == 8){

    set.seed(123)
    val <- alpha_nu(1, nu1 =2, 1, 2)

    expect_equal(round(val, 4), 0.6657)
  }
})

test_that("Expected value for alpha_nu when nu1 == 0",{
  set.seed(123)
  val <- alpha_nu(1, nu1 = 0, 1, 2)
  expect_equal(val, 0)
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

  expect_equal(prior.LST.R(c(1, 2), 1, 1, 2, T),
               prior_LST(c(1, 2), 1, 1, 2, T))
  expect_equal(prior.LST.R(c(1, 2), 1, 1, 2, F),
               prior_LST(c(1, 2), 1, 1, 2, F))
})

test_that("alpha_nu same in C++ as in R",{
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
                   (prior_nu(nu1, prior) / prior_nu(nu0, prior)))
    }
    return(aux)
  }

  expect_equal(alpha.nu(1, 2, c(1, 2), 2), alpha_nu(1, 2, c(1, 2), 2))

})

test_that("MH_nu_LST same in C++ as in R",{
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

  set.seed(123)
  MH.nu.LST.R <- MH.nu.LST(N = 1, omega2 = 0.5, beta = c(1, 2),
                           lambda = c(3, 4), nu0 = 0.2, prior = 2)

  set.seed(123)
  MH.nu.LST.Cpp <- MH_nu_LST(omega2 = 0.5, beta = c(1, 2),
                             lambda = c(3, 4), nu0 = 0.2, prior = 2)

  expect_equal(MH.nu.LST.R, MH.nu.LST.Cpp)

})

test_that("Post_lambda_obs_LST same in C++ as in R",{

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

  LST <- MCMC_LST(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])

  R.result <- Post.lambda.obs.LST(1, 1, cancer[, 3:11], LST)
  Cpp.result <- Post_lambda_obs_LST(1, 1, cancer[, 3:11], LST)
  expect_equal(R.result, Cpp.result)
})

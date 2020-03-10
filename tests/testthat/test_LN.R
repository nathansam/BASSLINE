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


test_that("MCMC_LN same in C++ as in R", {
  MCMC_LN_old <- function(N, thin, burn, Time, Cens, X, beta0 = NULL,
                          sigma20 = NULL, prior = 2, set = 1, eps_l = 0.5,
                          eps_r = 0.5) {

    # Sample starting values if not given
    if (is.null(beta0)) beta0 <- beta.sample(n = ncol(X))
    if (is.null(sigma20)) sigma20 <- sigma2.sample()

    MCMC.param.check(N, thin, burn, Time, Cens, X, beta0, sigma20,
                     prior, set, eps_l, eps_r)

    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- log(Time)

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    logt.aux <- logt[1, ]

    for (iter in 2:(N + 1)) {
        mu.aux <- solve(t(X) %*% X) %*% t(X) %*% logt.aux
        Sigma.aux <- sigma2.aux * solve(t(X) %*% X)
        beta.aux <- mvrnormArma(n = 1, mu = mu.aux, Sigma = Sigma.aux)


        shape.aux <- (n + 2 * p - 2) / 2
        rate.aux <- 0.5 * t(logt.aux - X %*% beta.aux) %*% (logt.aux - X %*%
                                                                beta.aux)
        if (rate.aux > 0 & is.na(rate.aux) == FALSE) {
            sigma2.aux <- (stats::rgamma(1, shape = shape.aux,
                                         rate = rate.aux)) ^ (-1)
        }


        logt.aux <- logt_update_SMLN(Time, Cens, X, beta.aux, sigma2.aux, set,
                                     eps_l, eps_r)

        if (iter %% thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            sigma2[iter / thin + 1] <- sigma2.aux
            logt[iter / thin + 1, ] <- logt.aux
        }
        if ((iter - 1) %% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    chain <- cbind(beta, sigma2, logt)



    beta.cols <- paste("beta.", seq(ncol(X)), sep = "")
    logt.cols <- paste("logt.", seq(length(Time)), sep = "")
    colnames(chain) <- c(beta.cols, "sigma2", logt.cols)

    if (burn > 0) {
      burn.period <- 1:(burn / thin)
      chain <- chain [-burn.period, ]
    }

    return(chain)
  }

#
#   set.seed(123)
#   Result.Cpp <- MCMC_LN(N = 1000, thin = 20, burn = 40,
#                         Time = cancer[, 1], Cens = cancer[, 2],
#                         X = cancer[, 3:11])
#
#   set.seed(123)
#   Result.R <- MCMC_LN_old(N = 1000, thin = 20, burn = 40,
#                                         Time = cancer[, 1], Cens = cancer[, 2],
#                                         X = cancer[, 3:11])
#
#   testthat::expect_equal(Result.Cpp, Result.R)



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
  Time <- seq(4)
  Cens <- c(0, 1, 0, 1)
  X <- matrix(seq(16), nrow = 4)
  beta <- c(2, 3, 4, 5)
  sigma2 <- 1.5; eps_l <- 0.5; eps_r <- 0.5

  for (set in 0:1){

    set.seed(123)
    result.R <- logt.update.SMLN(Time, Cens, X, beta, sigma2, set, eps_l,
                                 eps_r)
    set.seed(123)
    result.Cpp <- logt_update_SMLN(Time, Cens, X, beta, sigma2, set, eps_l,
                                   eps_r)

    expect_equal(as.numeric(result.R), as.numeric(result.Cpp))
  }
})



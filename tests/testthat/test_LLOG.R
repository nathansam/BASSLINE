################################################################################
################################ USER FUNCTIONS ################################
################################################################################
test_that("MCMC_LLOG returns expected num of rows when burn = 5 ,thin = 10",{

  N <- 100
  thin <- 10
  burn <- 5
  LLOG <- MCMC_LLOG(N = N, thin = thin, burn = burn, Time = cancer[,1],
                    Cens = cancer[,2], X = cancer[,3:11])
  expect_equal(nrow(LLOG), N / thin + 1 - burn)
})

test_that("RS_lambda_obs_LLOG returns same result in C++ as in R",{

  RS.lambda.obs.LLOG <- function(logt, X, beta, sigma2, obs, N.AKS) {
    lambda <- 0
    OK <- 0
    step <- 0
    while (OK == 0) {
      step <- step + 1

      lambda <- r_GIG(r = abs(logt[obs] - X[obs, ] %*% beta) / sqrt(sigma2))
      if (lambda != 0 && lambda != Inf) {
        U <- stats::runif(1, min = 0, max = 1)
        if (lambda > 4 / 3) {
          Z <- 1
          W <- exp(-0.5 * lambda)
          n.AKS <- 0
          while (n.AKS <= N.AKS) {
            n.AKS <- n.AKS + 1
            Z <- Z - ((n.AKS + 1) ^ 2) * W ^ (((n.AKS + 1) ^ 2) - 1)
            if (Z > U) {
              OK <- 1
            }
            n.AKS <- n.AKS + 1
            Z <- Z + ((n.AKS + 1) ^ 2) * W ^ (((n.AKS + 1) ^ 2) - 1)
            if (Z < U) {
              OK <- 0
            }
          }
        } else {
          H <- 0.5 * log(2) + 2.5 * log(pi) - 2.5 *
            log(lambda) - ((pi^2) / (2 * lambda)) + 0.5 *
            lambda
          lU <- log(U)
          Z <- 1
          W <- exp(- (pi ^ 2) / (2 * lambda))
          K <- lambda / (pi ^ 2)
          n.AKS <- 0
          while (n.AKS <= N.AKS) {
            n.AKS <- n.AKS + 1
            Z <- Z - K * W ^ ((n.AKS ^ 2) - 1)
            aux <- H + log(Z)
            if (is.na(aux) == FALSE && aux != -Inf && aux == Inf) {
              OK <- 1
            }
            if (is.na(aux) == FALSE && aux < Inf && aux > -Inf && aux > lU) {
              OK <- 1
            }
            n.AKS <- n.AKS + 1
            Z <- Z + ((n.AKS + 1) ^ 2) * W ^ (((n.AKS + 1) ^ 2) - 1)
            aux <- H + log(Z)
            if (is.na(aux) == FALSE && aux == -Inf) {
              OK <- 0
            }
            if (is.na(aux) == FALSE && aux > -Inf && aux < Inf && aux < lU) {
              OK <- 0
            }
          }
        }
      }
    }
    list(lambda = lambda, steps = step)
  }

  set.seed(123)
  R.result <- RS.lambda.obs.LLOG(c(1, 2), matrix(seq(4), ncol = 2), c(1, 2),
                                 0.2, 1, 1)

  set.seed(123)
  Cpp.result <- RS.lambda.obs.LLOG(c(1, 2), matrix(seq(4), ncol = 2), c(1, 2),
                                   0.2, 1, 1)

  expect_equal(R.result, Cpp.result)
})

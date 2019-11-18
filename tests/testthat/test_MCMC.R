

test_that("MCMC_LN first value",{
  data(cancer)
  set.seed(1)
  n <- dim(cancer)[1]
  Intercept <- rep(1, times = n) # Intercept
  x1 <- rep(0, times = n) # Treat (1 dummy variable)
  for(i in 1 : n) {if(cancer$treat[i] == 2) x1[i] <- 1}

  # Type (3 dummy variables)
  x2 <- rep(0, times = n)
  x3 <- rep(0, times = n)
  x4 <- rep(0, times = n)

  for(i in 1 : n) {if(cancer$type[i] == 1) x2[i] <- 1
  if(cancer$type[i] == 2) x3[i] <- 1
  if(cancer$type[i] == 3) x4[i] <- 1}
  x5 <- cancer$status # Status
  x6 <- cancer$mfd # Moths from diagnosis
  x7 <- cancer$age # Age
  x8 <- rep(0, times = n) # Prior (Dummy variable)
  for(i in 1 : n) {if(cancer$prior[i] == 10) x8[i] <- 1}
  X <- cbind(Intercept, x1, x2, x3, x4, x5, x6, x7, x8)
  Time <- cancer$time; Cens <- cancer$censor
  beta0 <- rnorm(9, 0, 1)
  sigma20 <- rgamma(1, 2, 2)

  LN <- MCMC_LN(N = 1000, thin = 20, Time, Cens, X, beta0, sigma20,
                prior = 2, set = 1, eps_l = 0.5, eps_r = 0.5)

  expect_equivalent(LN[1,1], -0.626453811)
})

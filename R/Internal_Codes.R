
################################################################################
#####################################PRIORS#####################################
################################################################################

###################################################### PRIOR FOR (BETA,SIGMA2)

#prior.LN <- function(beta, sigma2, prior, log) {
#    k = length(beta)
#    if (prior == 1) {
#        p <- 1 + k / 2
#    }
#    if (prior == 2) {
#        p <- 1
#    }
#    if (prior == 3) {
#        p <- 1
#   }
#    if (log == FALSE) {
#        aux <- sigma2 ^ (-p)
#    }
#    if (log == TRUE) {
#        aux <- -p * log(sigma2)
#    }
#    return(aux)
#}

#########################PRIOR FOR (BETA,SIGMA2,NU) (LOG-STUDENT'S T MODEL ONLY)

#IIJ.nu <- function(nu) {
#    aux <- sqrt(nu/ (nu + 3)) * sqrt(trigamma(nu / 2) -
#                                       trigamma((nu + 1) / 2) -
#                                            (2 * (nu + 3))/(nu * (nu + 1) ^ 2))
#    return(aux)
#}
prior.nu <- function(nu, prior) {
    if (prior == 2){
        aux <- IIJ_nu(nu) / stats::integrate(IIJ_nu, lower = 0,
                                             upper = Inf)$value
        return(aux)
    }
}
prior.LST <- function(beta, sigma2, nu, prior, log) {
    if (log == FALSE) {
        aux <- prior_LN(beta, sigma2, prior, logs = FALSE) *
          prior.nu(nu, prior)
    }
    if (log == TRUE) {
        aux <- prior_LN(beta, sigma2, prior, logs = TRUE) +
          log(prior.nu(nu, prior))
    }
    return(aux)
}

############### PRIOR FOR (BETA,SIGMA2,ALPHA) (LOG-EXPONENTIAL POWER MODEL ONLY)

#J.alpha <- function(alpha, k) {
#    aux <- ((alpha * (alpha - 1) * gamma(1 - 1 / alpha) /
#               gamma(1 / alpha)) ^ (k / 2)) * (1 / alpha) *
#                  sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
#    return(aux)
#}
#II.alpha <- function(alpha) {
#    aux <- (1 / alpha) * sqrt((1 + 1 / alpha) * trigamma(1 + 1 / alpha) - 1)
#    return(aux)
#}
#I.alpha <- function(alpha) {
#    aux <- sqrt(1 / alpha ^ 3) * sqrt((1 + 1 / alpha) *
#                 trigamma(1 + 1 / alpha) + (1 + digamma(1 + 1 / alpha))^ 2  - 1)
#    return(aux)
#}
prior.alpha <- function(alpha, k, prior) {
    if (prior == 1)
        aux <- J_alpha(alpha, k) / stats::integrate(J_alpha, lower = 1,
                                                    upper = 2, k = k)$value
    if (prior == 2)
        aux <- II_alpha(alpha) / stats::integrate(II_alpha, lower = 1,
                                                  upper = 2)$value
    if (prior == 3)
        aux <- I_alpha(alpha) / stats::integrate(I_alpha, lower = 1,
                                                 upper = 2)$value
    return(aux)
}
prior.LEP <- function(beta, sigma2, alpha, prior, log) {
    if (log == FALSE) {
        aux <- prior_LN(beta, sigma2, prior, logs = FALSE) *
                  prior.alpha(alpha, length(beta), prior)
    }
    if (log == TRUE) {
        aux <- prior_LN(beta, sigma2, prior, logs = TRUE) +
                  log(prior.alpha(alpha, length(beta), prior))
    }
    return(aux)
}


############# DATA AUGMENTATION UPDATE FOR RIGHT-CENSORED AND SET OBSERVATIONS
############# (LOGARITHMIC SCALE)

####################################################### FOR SMLN REPRESENTATION

logt.update.SMLN <- function(Time, Cens, X, beta, sigma2, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta
    if (set == 1) {
        aux <- Cens * (I(Time > eps_l) *
                         truncnorm::rtruncnorm(n = n,
                                               a = log(abs(Time - eps_l)),
                                               b = log(Time + eps_r),
                                               mean = MEAN,
                                               sd = sqrt(sigma2)) +
                  (1 - I(Time > eps_l)) *
                    truncnorm::rtruncnorm(n = n, a = -Inf,
                                          b = log(Time + eps_r),
                                          mean = MEAN,
                                          sd = sqrt(sigma2))) + (1 - Cens) *
                                    truncnorm::rtruncnorm(n = n,
                                                          a = log(Time),
                                                          b = Inf,
                                                          mean = MEAN,
                                                          sd = sqrt(sigma2))
    }
    if (set == 0) {
        aux <- Cens * log(Time) + (1 - Cens) *
                 truncnorm::rtruncnorm(n = n, a = log(Time), b = Inf,
                                       mean = MEAN, sd = sqrt(sigma2))
    }
    return(aux)
}

###### FOR MIXTURE OF UNIFORMS REPRESENTATION (LOG-EXPONENTIAL POWER MODEL ONLY)

logt.update.LEP <- function(Time, Cens, X, beta, sigma2, alpha, u, set, eps_l,
                            eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta

    if (set == 1) {
        a <- apply(cbind(MEAN - sqrt(sigma2) * u ^ (1 / alpha),
                         log(abs(Time - eps_l))), 1, max)
        a1 <- as.vector(MEAN - sqrt(sigma2) * u ^ (1 / alpha))
        b <- apply(cbind(MEAN + sqrt(sigma2) * u ^ (1 / alpha),
                         log(Time + eps_r)), 1, min)
        acens <- apply(cbind(MEAN - sqrt(sigma2) * u ^ (1 / alpha), log(Time)),
                       1, max)
        bcens <- as.vector(MEAN + sqrt(sigma2) * u ^ (1 / alpha))
        aux <- Cens * (as.numeric(I(Time > eps_l)) * (as.numeric(I(a < b)) *
                                      stats::runif(n, min = apply(cbind(a, b),
                                                                  1, min),
                                                   max = apply(cbind(a, b), 1,
                                                               max)) +
                                    (1 - as.numeric(I(a < b))) *
                                      log(Time)) +
                                       (1 - as.numeric(I(Time > eps_l))) *
                                        (as.numeric(I(a1 < b))
                                          * stats::runif(n,
                                                         min = apply(cbind(a1,
                                                                           b),
                                                                     1, min),
                                                         max = apply(cbind(a1,
                                                                           b),
                                                                     1, max)) +
                                           (1 - as.numeric(I(a1 < b))) *
                                           log(Time))) + (1 - Cens) *
                                           (as.numeric(I(acens < bcens)) *
            stats::runif(n, min = apply(cbind(acens, bcens), 1, min),
                         max = apply(cbind(acens, bcens), 1, max)) +
              (1 - as.numeric(I(acens < bcens))) * log(Time))
    }

    if (set == 0) {
        a <- apply(cbind(MEAN - sqrt(sigma2) * u^(1/alpha), log(Time)), 1, max)
        b <- as.vector(MEAN + sqrt(sigma2) * u^(1/alpha))
        aux <- Cens * log(Time) + (1 - Cens) * (as.numeric(I(a < b)) *
                      stats::runif(n, min = apply(cbind(a, b), 1, min),
                                   max = apply(cbind(a, b), 1, max)) +
                        (1 - as.numeric(I(a < b))) * log(Time))
    }
    return(aux)
}

###### OTHER FUNCTIONS REQUIRED FOR THE IMPLEMENTATION OF THE LOG-NORMAL MODEL
###### (NO MIXTURE)

######LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LN FUNCTIONS)

log.lik.LN <- function(Time, Cens, X, beta, sigma2, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, times = n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    if (set == 1) {
        # IF VERY EXTREME VALUES OF Time WERE OBSERVED,
        #  A NUMERICAL PROBLEM OCCURRED.  WE APPROXIMATED THE
        # LOG-LIKELIHOOD BY THE AREA OF A SQUARE UNDER THE CURVE.
        aux.set <- stats::plnorm(Time + eps_r, meanlog = MEAN,
                                 sdlog = sqrt(sigma2)) -
                      stats::plnorm(abs(Time - eps_l), meanlog = MEAN,
                                    sdlog = sqrt(sigma2))
        # ADD 2 TO AVOID -INF IN LOG (DOES NOT AFFECT THE RESULT)
        aux.set1 <- aux.set + I(aux.set == 0) * 2
        aux <- Cens * (I(Time > eps_l) * ((1 - I(aux.set == 0)) *
                                            log(aux.set1) + I(aux.set == 0) *
                                            (log(eps_l + eps_r) +
                                             stats::dlnorm(abs(Time - eps_l),
                                                           meanlog = MEAN,
                                                           sdlog = sqrt(sigma2),
                                                           log = TRUE))) +
                         (1 - I(Time > eps_l)) *
                         log(stats::plnorm(Time + eps_r, meanlog = MEAN,
                                           sdlog = sqrt(sigma2)) - 0)) +
                                                     (1 - Cens) *
                                       as.vector(stats::pnorm(MEAN - log(Time),
                                                              mean = 0,
                                                              sd = sqrt(sigma2),
                                                              log = TRUE))
    }
    if (set == 0) {
        aux <- Cens * stats::dlnorm(Time, meanlog = MEAN, sdlog = sqrt(sigma2),
                                    log = TRUE) + (1 - Cens) *
                       log(1 - stats::plnorm(Time, meanlog = MEAN,
                                             sdlog = sqrt(sigma2)))
    }
    return(sum(aux))
}

###########REDUCED CHAIN GIVEN A FIXED VALUE OF SIGMA2 (REQUIRED FOR ML.LN ONLY)

MCMCR.sigma2.LN <- function(N, thin, Time, Cens, X, beta0, sigma20, logt0,
                            prior, set, eps_l = 0.5, eps_r = 0.5) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k/2
    }
    if (prior == 2) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(sigma20, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    logt.aux <- logt[1, ]

    for (iter in 2:(N + 1)) {
        mu.aux <- solve(t(X) %*% X) %*% t(X) %*% logt.aux
        Sigma.aux <- sigma2.aux * solve(t(X) %*% X)
        beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux, sigma2.aux, set,
                                    eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            logt[iter/thin + 1, ] <- logt.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    chain <- cbind(beta, logt)
    return(chain)
}

### OTHER FUNCTIONS REQUIRED FOR THE IMPLEMENTATION OF THE LOG-STUDENT'S T MODEL

########METROPOLIS-HASTINMCMC UPDATE OF NU (REQUIRED FOR SEVERAL .LST FUNCTIONS)

MH.nu.LST <- function(N = 1, omega2, beta, lambda, nu0, prior) {
    k <- length(beta)
    ind <- rep(0, N)
    nu <- rep(0, times = N + 1)
    nu[1] <- nu0

    for (j.nu in 1:N) {
        y <- stats::rnorm(1, mean = nu[j.nu], sd = sqrt(omega2))
        ind.aux <- I(y >= 0)
        y <- abs(y)
        u.aux <- stats::runif(1, min = 0, max = 1)
        log.aux <- sum(stats::dgamma(lambda, shape = y/2, rate = y/2,
                                     log = TRUE)) -
                                  sum(stats::dgamma(lambda,
                                                    shape = nu[j.nu] / 2,
                                                    rate = nu[j.nu] / 2,
                                                    log = TRUE)) +
                                          log(prior.nu(y, prior) /
                                                prior.nu(nu[j.nu], prior))
        aux <- I(log(u.aux) < log.aux)
        aux <- as.numeric(aux) * as.numeric(ind.aux)
        nu[j.nu + 1] <- aux * y + (1 - aux) * nu[j.nu]
        ind[j.nu + 1] <- as.numeric(aux)
    }
    nu <- nu[-1]
    ind <- ind[-1]
    list(nu = nu, ind = ind)
}

#### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF NU
#### (REQUIRED FOR ML.LST FUNCTION ONLY)

alpha.nu <- function(nu0, nu1, lambda, k, prior) {
    if (nu1 <= 0) {
        aux <- 0
    }
    if (nu1 > 0) {
        aux <- min(1, exp(sum(stats::dgamma(lambda, shape = nu1 / 2,
                                            rate = nu1/2, log = TRUE))
                          - sum(stats::dgamma(lambda,
                                              shape = nu0 / 2,
                                              rate = nu0 / 2,
                                              log = TRUE))) *
                            (prior.nu(nu1, prior) / prior.nu(nu0, prior)))
    }
    return(aux)
}

################## LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LST FUNCTIONS)

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
                              log(1 - stats::pt((log(Time) - MEAN)/sqrt(sigma2),
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

#### NON-ADAPTIVE VERSION OF THE METROPOLIS-WITHIN-GIBSS ALGORITHM WITH FIXED
#### VALUES FOR THE VARIANCES OF THE GAUSSIAN PROPOSALS (REQUIRED FOR ML.LST
#### FUNCTION ONLY)

MCMC.LST.NonAdapt <- function(N, thin, Q, Time, Cens, X, beta0, sigma20, nu0,
                             prior, set, eps_l, eps_r, omega2.nu) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N/thin, 0)
    if (prior == 1) {
        p <- 1 + k/2
    }
    if (prior == 2) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    nu <- rep(0, times = N.aux + 1)
    nu[1] <- nu0
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- log(Time)
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- stats::rgamma(n, shape = nu0 / 2, rate = nu0/2)
    accept.nu <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    nu.aux <- nu[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]

    for (iter in 2:(N + 1)) {
        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        shape.aux <- (n + 2 * p - 2)/2
        rate.aux <- 0.5 * t(logt.aux - X %*% beta.aux) %*% Lambda %*%
                             (logt.aux - X %*% beta.aux)
        if (rate.aux > 0 & is.na(rate.aux) == FALSE) {
            sigma2.aux <- (stats::rgamma(1, shape = shape.aux,
                                         rate = rate.aux)) ^ (-1)
        }

        MH.nu <- MH.nu.LST(N = 1, omega2 = omega2.nu, beta.aux, lambda.aux,
                           nu.aux, prior)
        nu.aux <- MH.nu$nu
        accept.nu <- accept.nu + MH.nu$ind

        if ((iter - 1)%%Q == 0) {
            shape1.aux <- (nu.aux + 1)/2
            rate1.aux <- 0.5 * (nu.aux + ((logt.aux - X %*% beta.aux) ^ 2) /
                                  sigma2.aux)
            lambda.aux <- stats::rgamma(n, shape = rep(shape1.aux, times = n),
                                        rate = rate1.aux)
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                     sigma2.aux / lambda.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            sigma2[iter/thin + 1] <- sigma2.aux
            nu[iter/thin + 1] <- nu.aux
            logt[iter/thin + 1, ] <- logt.aux
            lambda[iter/thin + 1, ] <- lambda.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR nu :", round(accept.nu/N, 2), "\n"))

    chain <- cbind(beta, sigma2, nu, lambda, logt)
    return(chain)
}

############# REDUCED CHAIN GIVEN A FIXED VALUE OF NU (REQUIRED FOR ML.LST ONLY)

MCMCR.nu.LST <- function(N, thin, Q, Time, Cens, X, beta0, sigma20, nu0, logt0,
                         lambda0, prior, set, eps_l, eps_r) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k/2
    }
    if (prior == 2) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    nu <- rep(nu0, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- lambda0
    accept.nu <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    nu.aux <- nu[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]

    for (iter in 2:(N + 1)) {
        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        shape.aux <- (n + 2 * p - 2)/2
        rate.aux <- 0.5 * t(logt.aux - X %*% beta.aux) %*% Lambda %*%
                              (logt.aux - X %*% beta.aux)
        if (rate.aux > 0 & is.na(rate.aux) == FALSE) {
            sigma2.aux <- (stats::rgamma(1, shape = shape.aux,
                                         rate = rate.aux)) ^ (-1)
        }

        if ((iter - 1)%%Q == 0) {
            shape1.aux <- (nu.aux + 1)/2
            rate1.aux <- 0.5 * (nu.aux + ((logt.aux - X %*% beta.aux) ^ 2) /
                                  sigma2.aux)
            lambda.aux <- stats::rgamma(n, shape = rep(shape1.aux, times = n),
                                        rate = rate1.aux)
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                     sigma2.aux / lambda.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            sigma2[iter/thin + 1] <- sigma2.aux
            logt[iter/thin + 1, ] <- logt.aux
            lambda[iter/thin + 1, ] <- lambda.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR nu :", round(accept.nu / N, 2), "\n"))

    chain <- cbind(beta, sigma2, lambda, logt)
    return(chain)
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF SIGMA2 AND NU
#### (REQUIRED FOR ML.LST ONLY)

MCMCR.sigma2.nu.LST <- function(N, thin, Q, Time, Cens, X, beta0, sigma20, nu0,
                               logt0, lambda0, prior, set, eps_l, eps_r) {
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
    sigma2 <- rep(sigma20, times = N.aux + 1)
    nu <- rep(nu0, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- lambda0
    accept.nu <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    nu.aux <- nu[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]

    for (iter in 2:(N + 1)) {
        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        if ((iter - 1)%%Q == 0) {
            shape1.aux <- (nu.aux + 1)/2
            rate1.aux <- 0.5 * (nu.aux + ((logt.aux - X %*% beta.aux) ^ 2) /
                                  sigma2.aux)
            lambda.aux <- stats::rgamma(n, shape = rep(shape1.aux, times = n),
                                        rate = rate1.aux)
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                    sigma2.aux / lambda.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            logt[iter/thin + 1, ] <- logt.aux
            lambda[iter/thin + 1, ] <- lambda.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR nu :", round(accept.nu/N, 2), "\n"))

    chain <- cbind(beta, lambda, logt)
    return(chain)
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF LAMBDA[i]
#### (REQUIRED FOR BF.lambda.i.LST ONLY)

MCMCR.LST.lambda.obs <- function(ref, obs, N, thin, Q, Time, Cens, X, beta0,
                                 sigma20, nu0, lambda0, prior, set, eps_l,
                                 eps_r, ar = 0.44) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k/2
    }
    if (prior == 2) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    nu <- rep(0, times = N.aux + 1)
    nu[1] <- nu0
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- log(Time)
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- lambda0
    lambda[, obs] <- rep(ref, times = N.aux + 1)
    accept.nu <- 0
    pnu.aux <- 0
    ls.nu <- rep(0, times = N.aux + 1)

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    nu.aux <- nu[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]
    ls.nu.aux <- ls.nu[1]

    i_batch <- 0

    for (iter in 2:(N + 1)) {
        i_batch <- i_batch + 1

        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        shape.aux <- (n + 2 * p - 2)/2
        rate.aux <- 0.5 * t(logt.aux - X %*% beta.aux) %*% Lambda %*%
                               (logt.aux - X %*% beta.aux)
        if (rate.aux > 0 & is.na(rate.aux) == FALSE) {
            sigma2.aux <- (stats::rgamma(1, shape = shape.aux,
                                         rate = rate.aux)) ^ (-1)
        }

        MH.nu <- MH.nu.LST(N = 1, omega2 = exp(ls.nu.aux), beta.aux,
                           lambda.aux, nu.aux, prior)
        nu.aux <- MH.nu$nu
        accept.nu <- accept.nu + MH.nu$ind
        pnu.aux <- pnu.aux + MH.nu$ind

        if ((iter - 1)%%Q == 0) {
            shape1.aux <- (nu.aux + 1) / 2
            rate1.aux <- 0.5 * (nu.aux + ((logt.aux - X %*% beta.aux)^2) /
                                  sigma2.aux)
            lambda.aux <- stats::rgamma(n, shape = rep(shape1.aux, times = n),
                                        rate = rate1.aux)
            lambda.aux[obs] <- ref
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                     sigma2.aux / lambda.aux, set, eps_l, eps_r)


        if (i_batch == 50) {
            pnu.aux <- pnu.aux / 50
            Pnu.aux <- as.numeric(pnu.aux < ar)
            ls.nu.aux <- ls.nu.aux + ((-1)^Pnu.aux) * min(0.01, 1/sqrt(iter))
            i_batch <- 0
            pnu.aux <- 0
        }

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            sigma2[iter/thin + 1] <- sigma2.aux
            nu[iter/thin + 1] <- nu.aux
            logt[iter/thin + 1, ] <- logt.aux
            lambda[iter/thin + 1, ] <- lambda.aux
            ls.nu[iter/thin + 1] <- ls.nu.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR nu :", round(accept.nu/N, 2), "\n"))

    chain <- cbind(beta, sigma2, nu, lambda, logt, ls.nu)
    return(chain)
}

######## MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.LST ONLY)

Post.lambda.obs.LST <- function(obs, ref, X, chain) {
    N <- dim(chain)[1]
    n <- dim(X)[1]
    k <- dim(X)[2]
    aux1 <- rep(0, times = N)
    aux2 <- rep(0, times = N)

    for (iter in 1:N) {
        aux1[iter] <- (((chain[iter, (obs + k + 2 + n)] - X[obs, ] %*%
                           as.vector(chain[iter, (1:k)]))^2)/(chain[iter,
                                                                    (k + 1)])) +
                            chain[iter, (k + 2)]
        aux2[iter] <- stats::dgamma(x = ref,
                                    shape = (chain[iter, (k + 2)] + 1) / 2,
                                    rate = aux1[iter] / 2)
    }

    aux <- mean(aux2)
    return(aux)
}

#### CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF LAMBDA[obs]
#### (REQUIRED FOR BF.lambda.obs.LST ONLY)

CFP.obs.LST <- function(N, thin, Q, burn, ref, obs, Time, Cens, X, chain, prior,
                        set, eps_l, eps_r, ar = 0.44) {
    N.aux <- dim(chain)[1]
    k <- ncol(X)
    n <- nrow(X)
    chain1 <- MCMCR.LST.lambda.obs(ref, obs, N, thin, Q, Time, Cens, X,
                                   beta0 = t(chain[N.aux, 1:k]),
                                   sigma20 = chain[N.aux, (k + 1)],
                                   nu0 = chain[N.aux, (k + 2)],
                                   lambda0 = t(chain[N.aux,
                                                     (k + 3):(k + 2 + n)]),
                                   prior, set, eps_l, eps_r, ar)
    chain1 <- chain1[-(1:burn), ]
    N.aux2 <- dim(chain1)[1]
    aux1 <- rep(0, times = N.aux2)

    for (iter in 1:N.aux2) {
        aux1[iter] <- 1 / stats::dgamma(x = ref,
                                        shape = chain1[iter, (k + 2)] / 2,
                                        rate = chain1[iter, (k + 2)]/2)
    }
    aux <- mean(aux1)
    return(aux)
}

####### OTHER FUNCTIONS REQUIRED FOR THE IMPLEMENTATION OF THE LOG-LAPLACE MODEL

################# LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LLAP FUNCTIONS)

log.lik.LLAP <- function(Time, Cens, X, beta, sigma2, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, times = n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    if (set == 1) {
        aux <- Cens * (I(Time > eps_l) *
                         log(VGAM::plaplace(log(Time + eps_r),
                                            location = MEAN,
                                            scale = sqrt(sigma2)) -
                         VGAM::plaplace(log(abs(Time - eps_l)),
                                        location = MEAN,
                                        scale = sqrt(sigma2))) +
                           (1 - I(Time > eps_l)) *
                              log(VGAM::plaplace(log(Time + eps_r),
                                                 location = MEAN,
                                                 scale = sqrt(sigma2)) - 0)) +
          (1 - Cens) * log(1 - VGAM::plaplace(log(Time),
                                              location = MEAN,
                                              scale = sqrt(sigma2)))
    }
    if (set == 0) {
        aux <- Cens * (VGAM::dlaplace(log(Time), location = MEAN,
                                      scale = sqrt(sigma2), log = TRUE) -
                         log(Time)) + (1 - Cens) *
                           log(1 - VGAM::plaplace(log(Time),
                                                  location = MEAN,
                                                  scale = sqrt(sigma2)))
    }
    return(sum(aux))
}

######## REDUCED CHAIN GIVEN A FIXED VALUE OF SIGMA2 (REQUIRED FOR ML.LLAP ONLY)

MCMCR.sigma2.LLAP <- function(N, thin, Q, Time, Cens, X, beta0, sigma20, logt0,
                              lambda0, prior, set, eps_l = 0.5, eps_r = 0.5) {
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
    sigma2 <- rep(sigma20, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- lambda0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]

    for (iter in 2:(N + 1)) {
        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        if ((iter - 1)%%Q == 0) {
            mu.aux <- sqrt(sigma2.aux) / abs(logt.aux - X %*% beta.aux)
            if (sum(is.na(mu.aux)) == 0) {
                draw.aux <- VGAM::rinv.gaussian(n = n, mu = mu.aux,
                                                lambda = rep(1, times = n))
                lambda.aux <- I(draw.aux > 0) *
                                  draw.aux + (1 - I(draw.aux > 0)) *
                                      lambda.aux
            }
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                     sigma2.aux / lambda.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            logt[iter / thin + 1, ] <- logt.aux
            lambda[iter/thin + 1, ] <- lambda.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(iter)
        }
    }

    chain <- cbind(beta, lambda, logt)
    return(chain)
}

#### OTHER FUNCTIONS REQUIRED FOR THE IMPLEMENTATION OF THE
#### LOG-EXPONENTIAL POWER MODEL

#### METROPOLIS-HASTINMCMC UPDATE OF BETA[j]
#### (REQUIRED FOR SEVERAL .LEP FUNCTIONS)

MH.marginal.beta.j <- function(N = 1, omega2, logt, X, sigma2, alpha,
                               beta0, j) {
    k <- length(beta0)
    beta <- matrix(rep(0, times = k * (N + 1)), ncol = k)
    beta[1, ] <- beta0
    ind <- rep(0, times = N + 1)
    for (i.b in 1:N) {
        y.aux <- beta[i.b, ]
        y.aux[j] <- stats::rnorm(n = 1, mean = beta[i.b, j], sd = sqrt(omega2))
        aux <- -sum(abs((logt - X %*% y.aux) / sqrt(sigma2)) ^ alpha) +
                     sum(abs((logt - X %*% beta[i.b, ]) / sqrt(sigma2)) ^ alpha)
        u.aux <- stats::runif(1, 0, 1)
        if (is.na(aux)) {
            beta[i.b + 1, ] <- beta[i.b, ]
            ind[i.b + 1] <- 0
            cat("NA.beta\n")
        }
        if (is.na(aux) == FALSE && aux == -Inf) {
            beta[i.b + 1, ] <- beta[i.b, ]
            ind[i.b + 1] <- 0
            cat("-Inf.beta\n")
        }
        if (is.na(aux) == FALSE && aux == Inf) {
            beta[i.b + 1, ] <- y.aux
            ind[i.b + 1] <- 1
            cat("Inf.beta\n")
        }
        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
            beta[i.b + 1, ] <- y.aux
            ind[i.b + 1] <- 1
        }
        if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
            beta[i.b + 1, ] <- beta[i.b, ]
            ind[i.b + 1] <- 0
        }
    }
    beta <- beta[-1, ]
    ind <- ind[-1]
    list(beta = beta, ind = ind)
}

#### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF BETA[j]
#### (REQUIRED FOR ML.LEP FUNCTION ONLY)

alpha.beta <- function(beta.0, beta.1, logt, X, sigma2, alpha) {
    l1 <- -sum(abs((logt - X %*% beta.1) / sqrt(sigma2)) ^ alpha)
    l2 <- -sum(abs((logt - X %*% beta.0) / sqrt(sigma2)) ^ alpha)
    aux <- min(1, exp(l1 - l2))

    return(aux)
}

#### METROPOLIS-HASTINMCMC UPDATE OF SIGMA2
#### (REQUIRED FOR SEVERAL .LEP FUNCTIONS)

MH.marginal.sigma2 <- function(N = 1, omega2, logt, X, beta, alpha,
                               sigma20, prior) {
    k <- length(beta)
    n <- length(logt)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }
    sigma2 <- rep(0, times = N + 1)
    sigma2[1] <- sigma20
    ind <- rep(0, times = N + 1)
    for (i.s in 1:N) {
        y.aux <- stats::rnorm(n = 1, mean = sigma2[i.s], sd = sqrt(omega2))
        if (y.aux <= 0) {
            sigma2[i.s + 1] <- sigma2[i.s]
            ind[i.s + 1] <- 2
        } else {
            l1 <- -((n / 2) + p) * log(y.aux / sigma2[i.s])
            l2 <- -sum(abs(logt - X %*% beta) ^ alpha) * (y.aux ^ (-alpha / 2) -
                                                     sigma2[i.s] ^ (-alpha / 2))
            aux <- l1 + l2
            u.aux <- stats::runif(1, 0, 1)

            # THE FOLLOWING THREE LINES AVOID CRUSHES
            if (is.na(aux)) {
                sigma2[i.s + 1] <- sigma2[i.s]
                ind[i.s + 1] <- 0
                cat("NA.sigma2\n")
            }
            if (is.na(aux) == FALSE && aux == -Inf) {
                sigma2[i.s + 1] <- sigma2[i.s]
                ind[i.s + 1] <- 0
                cat("-Inf.sigma2\n")
            }
            if (is.na(aux) == FALSE && aux == Inf) {
                sigma2[i.s + 1] <- y.aux
                ind[i.s + 1] <- 1
                cat("Inf.sigma2\n")
            }

            if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
                sigma2[i.s + 1] <- y.aux
                ind[i.s + 1] <- 1
            }
            if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
                sigma2[i.s + 1] <- sigma2[i.s]
                ind[i.s + 1] <- 0
            }
        }
    }
    sigma2 <- sigma2[-1]
    ind <- ind[-1]
    list(sigma2 = sigma2, ind = ind)
}

#### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF SIGMA2
#### (REQUIRED FOR ML.LEP FUNCTION ONLY)

alpha.sigma2 <- function(sigma2.0, sigma2.1, logt, X, beta, alpha, prior) {
    if (sigma2.1 <= 0) {
        aux <- 0
    } else {
        k <- dim(X)[2]
        n <- dim(X)[1]
        if (prior == 1) {
            p <- 1 + k / 2
        }
        if (prior == 2) {
            p <- 1
        }
        if (prior == 3) {
            p <- 1
        }
        l1 <- -((n / 2) + p) * log(sigma2.1 / sigma2.0)
        l2 <- -sum(abs(logt - X %*% (beta))^alpha) *
                     (sigma2.1 ^ (-alpha / 2) - sigma2.0 ^ (-alpha / 2))
        aux <- min(1, exp(l1 + l2))
    }

    return(aux)
}

#### METROPOLIS-HASTINMCMC UPDATE OF ALPHA (REQUIRED FOR SEVERAL .LEP FUNCTIONS)

MH.marginal.alpha <- function(N = 1, omega2, logt, X, beta, sigma2,
                              alpha0, prior) {
    k <- dim(X)[2]
    n <- dim(X)[1]
    alpha <- rep(0, times = N + 1)
    alpha[1] <- alpha0
    ind <- rep(0, times = N + 1)

    for (i.a in 1:N) {
        y.aux <- stats::rnorm(n = 1, mean = alpha[i.a], sd = sqrt(omega2))
        if (y.aux >= 2 || y.aux <= 1) {
            alpha[i.a + 1] <- alpha[i.a]
            ind[i.a + 1] <- 0
        } else {
            l1 <- n * log(y.aux/gamma(1/y.aux)) -
                        sum(abs((logt - X %*% beta) / sqrt(sigma2)) ^ y.aux)
            l0 <- n * log(alpha[i.a] / gamma(1/alpha[i.a])) -
                       sum(abs((logt - X %*% beta) / sqrt(sigma2)) ^ alpha[i.a])
            aux <- l1 - l0 + log(prior.alpha(y.aux, k, prior) /
                                   prior.alpha(alpha[i.a], k, prior))
            u.aux <- stats::runif(1, 0, 1)

            # THE FOLLOWING THREE LINES AVOID CRUSHES
            if (is.na(aux)) {
                alpha[i.a + 1] <- alpha[i.a]
                ind[i.a + 1] <- 0
                cat("NA.alpha\n")
            }
            if (is.na(aux) == FALSE && aux == -Inf) {
                alpha[i.a + 1] <- alpha[i.a]
                ind[i.a + 1] <- 0
                cat("-Inf.alpha\n")
            }
            if (is.na(aux) == FALSE && aux == Inf) {
                alpha[i.a + 1] <- y.aux
                ind[i.a + 1] <- 1
                cat("Inf.alpha\n")
            }

            if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) < aux) {
                alpha[i.a + 1] <- y.aux
                ind[i.a + 1] <- 1
            }
            if (is.na(aux) == FALSE && aux > -Inf && log(u.aux) >= aux) {
                alpha[i.a + 1] <- alpha[i.a]
                ind[i.a + 1] <- 0
            }
        }
    }
    alpha <- alpha[-1]
    ind <- ind[-1]
    list(alpha = alpha, ind = ind)
}

#### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF ALPHA
#### (REQUIRED FOR ML.LEP FUNCTION ONLY)

alpha.alpha <- function(alpha0, alpha1, logt, X, beta, sigma2, prior) {
    k <- dim(X)[2]
    n <- dim(X)[1]
    if (alpha1 < 1 || alpha1 > 2) {
        aux <- 0
    } else {
        l1 <- n * log(alpha1 / gamma(1 / alpha1)) -
                 ((1 / sqrt(sigma2)) ^ alpha1) * sum((abs(logt - (X) %*%
                                                            (beta))) ^ alpha1)
        l0 <- n * log(alpha0 / gamma(1 / alpha0)) -
                 ((1 / sqrt(sigma2)) ^ alpha0) * sum((abs(logt - (X) %*%
                                                            (beta))) ^ alpha0)
        aux <- min(1, exp(l1 - l0) * prior.alpha(alpha1, k, prior) /
                                                 prior.alpha(alpha0, k, prior))
    }
    return(aux)
}

#### DISTRIBUTION FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
#### (BASED ON library normalp).

pnormp <- function(q, mu = 0, sigmap = 1, p = 2,
                   lower.tail = TRUE, log.pr = FALSE) {
 if (!is.numeric(q) || !is.numeric(mu) || !is.numeric(sigmap) || !is.numeric(p))
        stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1)
        stop("p must be at least equal to one")
    if (min(sigmap) <= 0)
        stop("sigmap must be positive")
    z <- (q - mu) / sigmap
    zz <- abs(z) ^ p
    zp <- stats::pgamma(zz, shape = 1 / p, scale = p)
    lzp <- stats::pgamma(zz, shape = 1 / p, scale = p, log = TRUE)
    zp <- ifelse(z < 0, 0.5 - exp(lzp - log(2)), 0.5 + exp(lzp - log(2)))
    if (log.pr == TRUE)
        zp <- ifelse(z < 0, log(0.5 - exp(lzp - log(2))),
                     log(0.5 + exp(lzp - log(2))))
    zp
}

#### DENSITY FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
#### (BASED ON library normalp).

dnormp <- function(x, mu = 0, sigmap = 1, p = 2, log = FALSE) {
 if (!is.numeric(x) || !is.numeric(mu) || !is.numeric(sigmap) || !is.numeric(p))
        stop(" Non-numeric argument to mathematical function")
    if (min(p) < 1)
        stop("p must be at least equal to one")
    if (min(sigmap) <= 0)
        stop("sigmap must be positive")
    cost <- 2 * p^(1 / p) * gamma(1 + 1 / p) * sigmap
    expon1 <- (abs(x - mu)) ^ p
    expon2 <- p * sigmap ^ p
    dsty <- (1 / cost) * exp(-expon1 / expon2)
    if (log == TRUE)
        dsty <- log(dsty)
    dsty
}

################## LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LEP FUNCTIONS)

log.lik.LEP <- function(Time, Cens, X, beta, sigma2, alpha, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    alpha <- rep(alpha, times = n)
    SP <- as.vector(sqrt(sigma2) * (1/alpha)^(1/alpha))
    if (set == 1) {
        aux1 <- (I(Time > eps_l) * log(pnormp(log(Time + eps_r),
                                              mu = MEAN,
                                              sigmap = SP,
                                              p = alpha) -
                                         pnormp(log(abs(Time - eps_l)),
                                                mu = MEAN, sigmap = SP,
                                                p = alpha)) +
                        (1 - I(Time > eps_l)) * log(pnormp(log(Time + eps_r),
                                                           mu = MEAN,
                                                           sigmap = SP,
                                                           p = alpha) - 0))
        aux2 <- log(1 - pnormp(log(Time), mu = MEAN, sigmap = SP, p = alpha))
        aux <- ifelse(Cens == 1, aux1, aux2)
    }
    if (set == 0) {
        aux <- Cens * (dnormp(log(Time), mu = MEAN, sigmap = SP,
                              p = alpha, log = TRUE) - log(Time)) + (1 - Cens) *
                     log(1 - pnormp(log(Time), mu = MEAN,
                                    sigmap = SP, p = alpha))
    }
    return(sum(aux))
}

#### NON-ADAPTIVE VERSION OF THE METROPOLIS-WITHIN-GIBSS ALGORITHM WITH
#### FIXED VALUES FOR THE VARIANCES OF THE GAUSSIAN PROPOSALS
#### (REQUIRED FOR ML.LEP FUNCTION ONLY)

MCMC.LEP.NonAdapt <- function(N, thin, Time, Cens, X, beta0, sigma20, alpha0,
                              prior, set, eps_l, eps_r, omega2.beta,
                              omega2.sigma2, omega2.alpha) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    alpha <- rep(0, times = N.aux + 1)
    alpha[1] <- alpha0
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- log(Time)
    U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    a <- ((abs(log(Time) - X %*% beta0)) / sqrt(sigma20))^alpha0
    U0 <- -log(1 - stats::runif(n)) + a
    U[1, ] <- U0

    accept.beta <- rep(0, times = k)
    accept.sigma2 <- 0
    accept.alpha <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    alpha.aux <- alpha[1]
    logt.aux <- logt[1, ]
    U.aux <- U[1, ]

    for (iter in 2:(N + 1)) {
        for (ind.b in 1:k) {
            MH.beta <- MH.marginal.beta.j(N = 1, omega2 = omega2.beta[ind.b],
                                          logt = logt.aux, X = X,
                                          sigma2 = sigma2.aux,
                                          alpha = alpha.aux,
                                          beta0 = beta.aux, j = ind.b)
            beta.aux[ind.b] <- MH.beta$beta[ind.b]
            if (MH.beta$ind == 1) {
                accept.beta[ind.b] <- accept.beta[ind.b] + 1
            }
        }

        MH.sigma2 <- MH.marginal.sigma2(N = 1, omega2 = omega2.sigma2,
                                        logt = logt.aux, X = X, beta = beta.aux,
                                        alpha = alpha.aux, sigma20 = sigma2.aux,
                                        prior = prior)
        sigma2.aux <- MH.sigma2$sigma2
        if (MH.sigma2$ind == 1) {
            accept.sigma2 <- accept.sigma2 + 1
        }

        MH.alpha <- MH.marginal.alpha(N = 1, omega2 = omega2.alpha,
                                      logt = logt.aux, X = X, beta = beta.aux,
                                      sigma2 = sigma2.aux, alpha0 = alpha.aux,
                                      prior = prior)
        alpha.aux <- MH.alpha$alpha
        if (MH.alpha$ind == 1) {
            accept.alpha <- accept.alpha + 1
        }

        a <- ((abs(logt.aux - X %*% beta.aux)) / sqrt(sigma2.aux)) ^ alpha.aux
        U.aux <- -log(1 - stats::runif(n)) + a

        logt.aux <- logt.update.LEP(Time, Cens, X, beta.aux, sigma2.aux,
                                    alpha.aux, u = U.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter/thin + 1, ] <- beta.aux
            sigma2[iter/thin + 1] <- sigma2.aux
            alpha[iter/thin + 1] <- alpha.aux
            logt[iter/thin + 1, ] <- logt.aux
            U[iter/thin + 1, ] <- U.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR beta", 1:k, ":", round(accept.beta/N, 2), "\n"))
    cat(paste("AR sigma2 :", round(accept.sigma2/N, 2), "\n"))
    cat(paste("AR alpha :", round(accept.alpha/N, 2), "\n"))

    chain = cbind(beta, sigma2, alpha, U, logt)
    return(chain)

}

########## REDUCED CHAIN GIVEN A FIXED VALUE OF ALPHA (REQUIRED FOR ML.LEP ONLY)

MCMCR.alpha.LEP <- function(N, thin, Time, Cens, X, beta0, sigma20, alpha0,
                            logt0, u0, prior, set, eps_l, eps_r,
                            omega2.beta, omega2.sigma2) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k/2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    alpha <- rep(alpha0, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    U[1, ] <- u0

    accept.beta <- rep(0, times = k)
    accept.sigma2 <- 0
    accept.alpha <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    alpha.aux <- alpha[1]
    logt.aux <- logt[1, ]
    U.aux <- U[1, ]

    for (iter in 2:(N + 1)) {
        for (ind.b in 1:k) {
            MH.beta <- MH.marginal.beta.j(N = 1, omega2 = omega2.beta[ind.b],
                                          logt = logt.aux, X = X,
                                          sigma2 = sigma2.aux,
                                          alpha = alpha.aux, beta0 = beta.aux,
                                          j = ind.b)
            beta.aux[ind.b] <- MH.beta$beta[ind.b]
            if (MH.beta$ind == 1) {
                accept.beta[ind.b] <- accept.beta[ind.b] + 1
            }
        }

        MH.sigma2 <- MH.marginal.sigma2(N = 1, omega2 = omega2.sigma2,
                                        logt = logt.aux, X = X, beta = beta.aux,
                                        alpha = alpha.aux, sigma20 = sigma2.aux,
                                        prior = prior)
        sigma2.aux <- MH.sigma2$sigma2
        if (MH.sigma2$ind == 1) {
            accept.sigma2 <- accept.sigma2 + 1
        }

        a <- ((abs(logt.aux - X %*% beta.aux)) / sqrt(sigma2.aux)) ^ alpha.aux
        U.aux <- -log(1 - stats::runif(n)) + a

        logt.aux <- logt.update.LEP(Time, Cens, X, beta.aux, sigma2.aux,
                                    alpha.aux, u = U.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            sigma2[iter / thin + 1] <- sigma2.aux
            logt[iter / thin + 1, ] <- logt.aux
            U[iter/thin + 1, ] <- U.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR beta", 1:k, ":", round(accept.beta/N, 2), "\n"))
    cat(paste("AR sigma2 :", round(accept.sigma2/N, 2), "\n"))
    cat(paste("AR alpha :", round(accept.alpha/N, 2), "\n"))

    chain <- cbind(beta, sigma2, U, logt)
    return(chain)
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF SIGMA2 AND ALPHA
#### (REQUIRED FOR ML.LEP ONLY)

MCMCR.sigma2.alpha.LEP <- function(N, thin, Time, Cens, X, beta0, sigma20,
                                   alpha0, logt0, u0, prior, set, eps_l,
                                   eps_r, omega2.beta) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(sigma20, times = N.aux + 1)
    alpha <- rep(alpha0, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    U[1, ] <- u0

    accept.beta <- rep(0, times = k)
    accept.sigma2 <- 0
    accept.alpha <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    alpha.aux <- alpha[1]
    logt.aux <- logt[1, ]
    U.aux <- U[1, ]

    for (iter in 2:(N + 1)) {
        for (ind.b in 1:k) {
            MH.beta <- MH.marginal.beta.j(N = 1, omega2 = omega2.beta[ind.b],
                                          logt = logt.aux, X = X,
                                          sigma2 = sigma2.aux,
                                          alpha = alpha.aux, beta0 = beta.aux,
                                          j = ind.b)
            beta.aux[ind.b] <- MH.beta$beta[ind.b]
            if (MH.beta$ind == 1) {
                accept.beta[ind.b] <- accept.beta[ind.b] + 1
            }
        }

        a <- ((abs(logt.aux - X %*% beta.aux)) / sqrt(sigma2.aux)) ^ alpha.aux
        U.aux <- -log(1 - stats::runif(n)) + a

        logt.aux <- logt.update.LEP(Time, Cens, X, beta.aux, sigma2.aux,
                                    alpha.aux, u = U.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            logt[iter / thin + 1, ] <- logt.aux
            U[iter / thin + 1, ] <- U.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
           cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR beta", 1:k, ":", round(accept.beta/N, 2), "\n"))
    cat(paste("AR sigma2 :", round(accept.sigma2/N, 2), "\n"))
    cat(paste("AR alpha :", round(accept.alpha/N, 2), "\n"))

    chain <- cbind(beta, U, logt)
    return(chain)
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF SOME BETA'S, SIGMA2 AND ALPHA
#### (REQUIRED FOR ML.LEP ONLY)

MCMCR.betaJ.sigma2.alpha.LEP <- function(N, thin, Time, Cens, X, beta0, sigma20,
                                         alpha0, logt0, u0, prior, set,
                                         eps_l, eps_r, omega2.beta, J) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(sigma20, times = N.aux + 1)
    alpha <- rep(alpha0, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    U[1, ] <- u0

    accept.beta <- rep(0, times = k)
    accept.sigma2 <- 0
    accept.alpha <- 0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    alpha.aux <- alpha[1]
    logt.aux <- logt[1, ]
    U.aux <- U[1, ]

    for (iter in 2:(N + 1)) {
        if (J < k) {
            for (ind.b in (J + 1):k) {
                MH.beta <- MH.marginal.beta.j(N = 1,
                                              omega2 = omega2.beta[ind.b],
                                              logt = logt.aux, X = X,
                                              sigma2 = sigma2.aux,
                                              alpha = alpha.aux,
                                              beta0 = beta.aux,
                                              j = ind.b)
                beta.aux[ind.b] <- MH.beta$beta[ind.b]
                if (MH.beta$ind == 1) {
                  accept.beta[ind.b] <- accept.beta[ind.b] + 1
                }
            }
        }

        a <- ((abs(logt.aux - X %*% beta.aux)) / sqrt(sigma2.aux)) ^ alpha.aux
        U.aux <- -log(1 - stats::runif(n)) + a

        logt.aux <- logt.update.LEP(Time, Cens, X, beta.aux, sigma2.aux,
                                    alpha.aux, u = U.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            logt[iter / thin + 1, ] <- logt.aux
            U[iter / thin + 1, ] <- U.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR beta", 1:k, ":", round(accept.beta / N, 2), "\n"))
    cat(paste("AR sigma2 :", round(accept.sigma2 / N, 2), "\n"))
    cat(paste("AR alpha :", round(accept.alpha / N, 2), "\n"))

    chain <- cbind(beta, U, logt)
    return(chain)
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF U[i]
#### (REQUIRED FOR BF.lambda.i.LEP ONLY)

MCMCR.LEP.u.i <- function(ref, obs, N, thin, Time, Cens, X, beta0, sigma20,
                         alpha0, u0, prior, set, eps_l, eps_r, ar = 0.44) {
    k <- length(beta0)
    n <- length(Time)
    N.aux <- round(N / thin, 0)
    if (prior == 1) {
        p <- 1 + k / 2
    }
    if (prior == 2) {
        p <- 1
    }
    if (prior == 3) {
        p <- 1
    }

    beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    beta[1, ] <- beta0
    sigma2 <- rep(0, times = N.aux + 1)
    sigma2[1] <- sigma20
    alpha <- rep(0, times = N.aux + 1)
    alpha[1] <- alpha0
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- log(Time)
    U <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    U[1, ] <- u0

    accept.beta <- rep(0, times = k)
    pbeta.aux <- rep(0, times = k)
    ls.beta <- matrix(rep(0, times = (N.aux + 1) * k), ncol = k)
    accept.sigma2 <- 0
    psigma2.aux <- 0
    ls.sigma2 <- rep(0, times = N.aux + 1)
    accept.alpha <- 0
    palpha.aux <- 0
    ls.alpha <- rep(0, times = N.aux + 1)

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    alpha.aux <- alpha[1]
    logt.aux <- logt[1, ]
    U.aux <- U[1, ]
    ls.beta.aux <- ls.beta[1, ]
    ls.sigma2.aux <- ls.sigma2[1]
    ls.alpha.aux <- ls.alpha[1]

    i_batch <- 0

    for (iter in 2:(N + 1)) {
        i_batch <- i_batch + 1

        for (ind.b in 1:k) {
            MH.beta <- MH.marginal.beta.j(N = 1,
                                          omega2 = exp(ls.beta.aux[ind.b]),
                                          logt = logt.aux, X = X,
                                          sigma2 = sigma2.aux,
                                          alpha = alpha.aux,
                                          beta0 = beta.aux,
                                          j = ind.b)
            beta.aux[ind.b] <- MH.beta$beta[ind.b]
            if (MH.beta$ind == 1) {
                accept.beta[ind.b] <- accept.beta[ind.b] + 1
                pbeta.aux[ind.b] <- pbeta.aux[ind.b] + 1
            }
        }

        MH.sigma2 <- MH.marginal.sigma2(N = 1, omega2 = exp(ls.sigma2.aux),
                                        logt = logt.aux, X = X, beta = beta.aux,
                                        alpha = alpha.aux, sigma20 = sigma2.aux,
                                        prior = prior)
        sigma2.aux <- MH.sigma2$sigma2
        if (MH.sigma2$ind == 1) {
            accept.sigma2 <- accept.sigma2 + 1
            psigma2.aux <- psigma2.aux + 1
        }

        MH.alpha <- MH.marginal.alpha(N = 1, omega2 = exp(ls.alpha.aux),
                                      logt = logt.aux, X = X, beta = beta.aux,
                                      sigma2 = sigma2.aux, alpha0 = alpha.aux,
                                      prior = prior)
        alpha.aux <- MH.alpha$alpha
        if (MH.alpha$ind == 1) {
            accept.alpha <- accept.alpha + 1
            palpha.aux <- palpha.aux + 1
        }

        a <- ((abs(logt.aux - X %*% beta.aux)) /
                sqrt(sigma2.aux)) ^ alpha.aux
        U.aux <- -log(1 - stats::runif(n)) + a
        U.aux[obs] <- ref

        logt.aux <- logt.update.LEP(Time, Cens, X, beta.aux, sigma2.aux,
                                    alpha.aux, u = U.aux, set, eps_l, eps_r)

        if (i_batch == 50) {
            pbeta.aux <- pbeta.aux / 50
            Pbeta.aux <- as.numeric(pbeta.aux < rep(ar, times = k))
            ls.beta.aux <- ls.beta.aux + ((-1) ^ Pbeta.aux) *
                                                    min(0.01, 1 / sqrt(iter))
            psigma2.aux <- psigma2.aux / 50
            Psigma2.aux <- as.numeric(psigma2.aux < ar)
            ls.sigma2.aux <- ls.sigma2.aux + ((-1) ^ Psigma2.aux) *
                                                    min(0.01, 1 / sqrt(iter))
            palpha.aux <- palpha.aux / 50
            Palpha.aux <- as.numeric(palpha.aux < ar)
            ls.alpha.aux <- ls.alpha.aux + ((-1) ^ Palpha.aux) *
                                                     min(0.01, 1 / sqrt(iter))
            i_batch <- 0
            pbeta.aux <- rep(0, times = k)
            psigma2.aux <- 0
            palpha.aux <- 0
        }

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            sigma2[iter / thin + 1] <- sigma2.aux
            alpha[iter / thin + 1] <- alpha.aux
            logt[iter / thin + 1, ] <- logt.aux
            U[iter / thin + 1, ] <- U.aux
            ls.beta[iter / thin + 1, ] <- ls.beta.aux
            ls.sigma2[iter / thin + 1] <- ls.sigma2.aux
            ls.alpha[iter / thin + 1] <- ls.alpha.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    cat(paste("AR beta", 1:k, ":", round(accept.beta / N, 2), "\n"))
    cat(paste("AR sigma2 :", round(accept.sigma2 / N, 2), "\n"))
    cat(paste("AR alpha :", round(accept.alpha / N, 2), "\n"))

    chain <- cbind(beta, sigma2, alpha, U, logt, ls.beta, ls.sigma2, ls.alpha)
    return(chain)

}

#### DENSITY FUNCTION OF A TRUNCATED EXPONENTIAL DISTRIBUTION
#### (REQUIRED FOR BF.u.i.LEP ONLY)

dtexp <- function(x, rate, trunc) {
    if (x >= trunc) {
        aux <- rate * exp(-rate * (x - trunc))
    } else {
        aux <- 0
    }
    return(aux)
}

#### MARGINAL POSTERIOR OF u[obs]
#### (REQUIRED FOR BF.u.obs.LEP ONLY)

Post.u.obs.LEP <- function(obs, ref, X, chain) {
    N <- dim(chain)[1]
    n <- dim(X)[1]
    k <- dim(X)[2]
    aux1 <- rep(0, times = N)
    aux2 <- rep(0, times = N)

    for (iter in 1:N) {
        trunc.aux <- (abs(chain[iter, (obs + k + 2 + n)] - X[obs, ] %*%
                            as.vector(chain[iter, 1:k])) /
                        sqrt(chain[iter, k + 1])) ^ (chain[iter, k + 2])

        aux1[iter] <- dtexp(x = ref, rate = 1, trunc = trunc.aux)
    }
    aux <- mean(aux1)
    return(aux)
}

#### CORRECTION FACTOR/PRIOR FOR BAYES FACTOR OF u[obs]
#### (REQUIRED FOR BF.u.obs.LEP ONLY)

CFP.obs.LEP <- function(N, thin, burn, ref, obs, Time, Cens, X, chain, prior,
                        set, eps_l, eps_r, ar = 0.44) {
    N.aux <- dim(chain)[1]
    k <- ncol(X)
    n <- nrow(X)
    chain1 <- MCMCR.LEP.u.i(ref, obs, N, thin, Time, Cens, X,
                            beta0 = as.vector(chain[N.aux, 1:k]),
                            sigma20 = chain[N.aux, (k + 1)],
                            alpha0 = chain[N.aux, (k + 2)],
                            u0 = as.vector(chain[N.aux, (k + 3):(k + 2 + n)]),
                            prior, set, eps_l, eps_r, ar)
    chain1 <- chain1[-(1:burn), ]
    N.aux2 <- dim(chain1)[1]
    aux1 <- rep(0, times = N.aux2)

    for (iter in 1:N.aux2) {
        aux1[iter] <- 1 / stats::dgamma(x = ref,
                                        shape = 1 + (1/chain1[iter, k + 2]),
                                        rate = 1)
    }

    aux <- mean(aux1)
    return(aux)
}

#### OTHER FUNCTIONS REQUIRED FOR THE IMPLEMENTATION OF THE LOG-LOGISTIC MODEL

#### SAMPLING FROM A GENERALIZED INVERSE GAUSSIAN DISTRIBUTION
#### (REQUIRED FOR SEVERAL .LLOG FUNCTIONS)
#### BASED ON HOLMES AND HELD (2006)

rGIG <- function(n = 1, r) {
    y <- stats::rnorm(n)
    y <- y^2
    y <- 1 + ((y - sqrt(y * (4 * r + y))) / (2 * r))
    u <- stats::runif(n)
    aux <- rep(0, times = n)
    for (i in 1:n) {
        if (u[i] <= 1/(1 + y[i])) {
            aux[i] <- r/y[i]
        } else {
            aux[i] <- r * y[i]
        }
    }
    return(aux)
}

#### REJECTION SAMPLING FOR LAMBDA[i]
#### (REQUIRED FOR SEVERAL .LLOG FUNCTIONS)
#### BASED ON HOLMES AND HELD (2006)

RS.lambda.obs.LLOG <- function(logt, X, beta, sigma2, obs, N.AKS) {
    lambda <- 0
    OK <- 0
    step <- 0
    while (OK == 0) {
        step <- step + 1

        lambda <- rGIG(n = 1,
                       r = abs(logt[obs] - X[obs, ] %*% beta) / sqrt(sigma2))
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
                W <- exp(-(pi ^ 2) / (2 * lambda))
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

#### LOG-LIKELIHOOD FUNCTION
#### (REQUIRED FOR SEVERAL .LLOG FUNCTIONS)

log.lik.LLOG <- function(Time, Cens, X, beta, sigma2, set, eps_l, eps_r) {
    n <- length(Time)
    aux <- rep(0, n)
    MEAN <- X %*% beta
    sigma2 <- rep(sigma2, times = n)
    if (set == 1) {
        aux <- Cens * (I(Time > eps_l) *
                         log(stats::plogis(log(Time + eps_r), location = MEAN,
                                           scale = sqrt(sigma2)) -
                               stats::plogis(log(abs(Time - eps_l)),
                                             location = MEAN,
                                             scale = sqrt(sigma2))) +
                                     (1 - I(Time > eps_l)) *
                                 log(stats::plogis(log(Time + eps_r),
                                                   location = MEAN,
                                                   scale = sqrt(sigma2)) - 0)) +
                                        (1 - Cens) *
                                 log(1 - stats::plogis(log(Time),
                                                       location = MEAN,
                                                       scale = sqrt(sigma2)))
    }
    if (set == 0) {
        aux <- Cens * (stats::dlogis(log(Time),
                                     location = MEAN,
                                     scale = sqrt(sigma2),
                                     log = TRUE) - log(Time)) + (1 - Cens) *
                          log(1 - stats::plogis(log(Time), location = MEAN,
                                                scale = sqrt(sigma2)))
    }
    return(sum(aux))
}

#### REDUCED CHAIN GIVEN A FIXED VALUE OF SIGMA2
#### (REQUIRED FOR ML.LLOG ONLY)

MCMCR.sigma2.LLOG <- function(N, thin, Q, Time, Cens, X, beta0, sigma20, logt0,
                              lambda0, prior, set, eps_l = 0.5, eps_r = 0.5,
                              N.AKS = 3) {
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
    sigma2 <- rep(sigma20, times = N.aux + 1)
    logt <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    logt[1, ] <- logt0
    lambda <- matrix(rep(0, times = (N.aux + 1) * n), ncol = n)
    lambda[1, ] <- lambda0

    beta.aux <- beta[1, ]
    sigma2.aux <- sigma2[1]
    logt.aux <- logt[1, ]
    lambda.aux <- lambda[1, ]

    for (iter in 2:(N + 1)) {
        Lambda <- diag(lambda.aux)
        AUX1 <- (t(X) %*% Lambda %*% X)
        if (det(AUX1) != 0) {
            AUX <- solve(AUX1)
            mu.aux <- AUX %*% t(X) %*% Lambda %*% logt.aux
            Sigma.aux <- sigma2.aux * AUX
            beta.aux <- MASS::mvrnorm(n = 1, mu = mu.aux, Sigma = Sigma.aux)
        }

        if ((iter - 1)%%Q == 0) {
            for (obs in 1:n) {
                lambda.aux[obs] <- 1 / RS.lambda.obs.LLOG(logt = logt.aux,
                                                          X = X,
                                                          beta = beta.aux,
                                                          sigma2 = sigma2.aux,
                                                          obs = obs,
                                                          N.AKS = N.AKS)$lambda
            }
        }

        logt.aux <- logt.update.SMLN(Time, Cens, X, beta.aux,
                                     sigma2.aux / lambda.aux, set, eps_l, eps_r)

        if (iter%%thin == 0) {
            beta[iter / thin + 1, ] <- beta.aux
            sigma2[iter / thin + 1] <- sigma2.aux
            logt[iter / thin + 1, ] <- logt.aux
            lambda[iter / thin + 1, ] <- lambda.aux
        }
        if ((iter - 1) %*% 1e+05 == 0) {
            cat(paste("Iteration :", iter, "\n"))
        }
    }

    chain <- cbind(beta, lambda, logt)
    return(chain)
}


#### Input check
MCMC.param.check <- function(N, thin, burn, Time, Cens, X, beta0, sigma20,
                             prior, set, eps_l, eps_r){

    num.obs <- nrow(X)
    num.covariates <- ncol(X)

    if (is.matrix(X) == F){
        stop("X is not a matrix.\n")
    }

    # Check N is a  0 < integer
    if (N <= 0 | N %% 1 != 0 ){
        stop("N should be an integer greater than zero.\n")
    }

    if (thin < 2 | thin %% 1 != 0){
        stop("thin should be a integer > 2.\n")
    }

    if (burn < 0 | burn %% 1 != 0){
        stop("burn should be a non-negative integer.\n")
    }

    if (N / thin  < burn ){
        stop("The number of non-discarded iterations will be zero!\n Adjust N, thin or burn.\n")
    }

    if (all (Time > 0) == F){
        stop ("All values in Time should be non-negative.\n")
    }
    if (length(Time) != num.obs){
    stop ("Time is not the correct length.\n")
        }

    if ( all(Cens == 1 | Cens == 0) == F){
        stop("Cens should be either 0 or 1 for each observation\n")
    }
    if (length(Cens) != num.obs){
        stop("Cens is not the correct length.\n")
    }

    if (length(beta0) != num.covariates){
        stop("beta0 is not the correct length.\n")
    }

    if (prior %in% c(1,2,3) == F){
        stop("prior should be 1, 2 or 3. See documentation\n")
    }

    if (set %in% c(1,2) == F){
        stop("set should be 1 or 2. See documentation\n")
    }

    if (burn == 0){cat(paste0("Note! No burn-in period is being used!\n"))}

}

beta.sample <- function(n){
    cat("Sampling initial betas from a Normal(0, 1) distribution\n")
    betas <- stats::rnorm(n, 0, 1)
    cat(paste("Initial beta", 1:length(betas), ":", round(betas, 2), "\n"))
    cat("\n")
    return(betas)
}

sigma2.sample <- function(){
    cat("Sampling initial sigma^2 from a Gamma(2, 2) distribution\n")
    sigma2 <- stats::rgamma(1, 2, 2)
    cat(paste("Initial sigma^2 :", round(sigma2, 2), "\n\n"))
    return(sigma2)
}

nu.sample <- function(){
    cat("Sampling initial nu from a Gamma(2, 2) distribution\n")
    nu <- stats::rgamma(1, 2, 2)
    cat(paste("Initial nu :", round(nu, 2), "\n\n"))
    return(nu)

}

alpha.sample <- function(){
    cat("Sampling initial alpha from a Uniform(1, 2) distribution\n")
    alpha <- stats::runif(1, 1, 2)
    cat(paste("Initial alpha :", round(alpha, 2), "\n\n"))
    return(alpha)
}

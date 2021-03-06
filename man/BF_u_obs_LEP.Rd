% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LogExponentialPower.R
\name{BF_u_obs_LEP}
\alias{BF_u_obs_LEP}
\title{Outlier detection for observation for the log-exponential power model}
\usage{
BF_u_obs_LEP(
  N,
  thin,
  burn,
  ref,
  obs,
  Time,
  Cens,
  X,
  chain,
  prior = 2,
  set = TRUE,
  eps_l = 0.5,
  eps_r = 0.5,
  ar = 0.44
)
}
\arguments{
\item{N}{Total number of iterations. Must be a multiple of thin.}

\item{thin}{Thinning period.}

\item{burn}{Burn-in period}

\item{ref}{Reference value \eqn{u_{ref}}. Vallejos & Steel recommends this
value be set to \eqn{1.6 +1_\alpha} for the LEP model.}

\item{obs}{Indicates the number of the observation under analysis}

\item{Time}{Vector containing the survival times.}

\item{Cens}{Censoring indication (1: observed, 0: right-censored).}

\item{X}{Design matrix with dimensions \eqn{n} x  \eqn{k} where \eqn{n} is
the number of observations and \eqn{k} is the number of covariates
(including the intercept).}

\item{chain}{MCMC chains generated by a BASSLINE MCMC function}

\item{prior}{Indicator of prior (1: Jeffreys, 2: Type I Ind. Jeffreys,
3: Ind. Jeffreys).}

\item{set}{Indicator for the use of set observations (1: set observations,
0: point observations). The former is strongly recommended over the
latter as point observations cause problems in the context of Bayesian
inference (due to continuous sampling models assigning zero probability
to a point).}

\item{eps_l}{Lower imprecision \eqn{(\epsilon_l)} for set observations
(default value: 0.5).}

\item{eps_r}{Upper imprecision \eqn{(\epsilon_r)} for set observations
(default value: 0.5)}

\item{ar}{Optimal acceptance rate for the adaptive Metropolis-Hastings
updates}
}
\description{
This returns a unique number corresponding to the Bayes Factor
associated to the test \eqn{M_0: \Lambda_{obs} = \lambda_{ref}} versus
\eqn{M_1: \Lambda_{obs}\neq \lambda_{ref}} (with all other
\eqn{\Lambda_j,\neq obs} free). The value of \eqn{\lambda_{ref}} is
required as input. The user should expect long running times for the
log-Student’s t model, in which case a reduced chain given
\eqn{\Lambda_{obs} = \lambda_{ref}} needs to be generated
}
\examples{
library(BASSLINE)

# Please note: N=1000 is not enough to reach convergence.
# This is only an illustration. Run longer chains for more accurate
# estimations (especially for the log-exponential power model).

LEP <- MCMC_LEP(N = 1000, thin = 20, burn = 40, Time = cancer[, 1],
                Cens = cancer[, 2], X = cancer[, 3:11])
alpha <- mean(LEP[, 11])
uref <- 1.6 + 1 / alpha
LEP.Outlier <- BF_u_obs_LEP(N = 100, thin = 20, burn =1 , ref = uref,
                            obs = 1, Time = cancer[, 1], Cens = cancer[, 2],
                            cancer[, 3:11], chain = LEP)

}

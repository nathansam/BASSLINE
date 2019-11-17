% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/User_Codes.R
\name{MCMC_LLOG}
\alias{MCMC_LLOG}
\title{MCMC algorithm for the log-logistic model}
\usage{
MCMC_LLOG(
  N,
  thin,
  Q,
  Time,
  Cens,
  X,
  beta0,
  sigma20,
  prior,
  set,
  eps_l = 0.5,
  eps_r = 0.5,
  N.AKS = 3
)
}
\arguments{
\item{N}{Total number of iterations.}

\item{thin}{Thinning period.}

\item{Q}{Update period for the \eqn{\lambda_{i}}’s}

\item{Time}{Vector containing the survival times.}

\item{Cens}{Censoring indication (1: observed, 0: right-censored).}

\item{X}{Design matrix with dimensions \eqn{n} x  \eqn{k} where \eqn{n} is
the number of observations and \eqn{k} is the number of covariates
(including the intercept).}

\item{beta0}{Starting value for \eqn{\beta}.}

\item{sigma20}{Starting value for \eqn{\sigma^2}.}

\item{prior}{Indicator of prior (1: Jeffreys, 2: Type I Ind. Jeffreys,
3: Ind. Jeffreys).}

\item{set}{Indicator for the use of set observations (1: set observations,
0: point observations).}

\item{eps_l}{Lower imprecision \eqn{(\epsilon_l)} for set observations
(default value: 0.5).}

\item{eps_r}{Upper imprecision \eqn{(\epsilon_r)} for set observations
(default value: 0.5)}

\item{N.AKS}{Maximum number of terms of the Kolmogorov-Smirnov density used
for the rejection sampling when updating mixing parameters (default
value: 3)}
}
\description{
Markov chain Monte carlo algorithm for the log-logistic model
}
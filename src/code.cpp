#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <RcppArmadillo.h>


using Rcpp::NumericVector;

/*##############################################################################
#####################################PRIORS#####################################
##############################################################################*/

//########################### PRIOR FOR (BETA,SIGMA2) ##########################

// [[Rcpp::export]]
double prior_LN(NumericVector beta, double sigma2, int prior, bool logs){

  float p;
  double aux;

  int k = beta.size();

  if (prior == 1){
    p = 1.0 + k / 2.0;
  }

  if (prior == 2){
    p = 1;
  }

  if (prior == 3){
    p = 1;
  }

  if (logs == false){
    aux = sigma2 * (-p);
  }

  if (logs == true){
    aux = log(sigma2) * -p;
  }
  return aux;
}


// [[Rcpp::export]]
NumericVector prior_nu(NumericVector nu, int prior){

  if (prior == 2){
    NumericVector aux = sqrt(nu / (nu + 3)) *
      sqrt (Rcpp::trigamma(nu / 2) - Rcpp::trigamma((nu + 1) / 2) -
      (2 * (nu + 3)) / (nu * pow(nu + 1,2)));
    return aux;
  }
  return NULL;
}


//########### PRIOR FOR (BETA,SIGMA2,NU) (LOG-STUDENT'S T MODEL ONLY)###########

// [[Rcpp::export]]
NumericVector prior_LST(NumericVector beta, double sigma2, NumericVector nu,
                        int prior, bool logs) {

  if (logs == false) {
    NumericVector aux = prior_LN(beta, sigma2, prior,
                                 false) * prior_nu(nu, prior);
    return aux;
  } else {
    NumericVector aux = prior_LN(beta, sigma2, prior,
                                 true) + log(prior_nu(nu, prior));
    return aux;
  }
}




// [[Rcpp::export]]
NumericVector J_alpha(NumericVector alpha, int k){
  NumericVector aux = pow(alpha * (alpha - 1) *

    Rcpp::trigamma(1.0 - 1.0 / alpha) /
      Rcpp::trigamma (1.0 / alpha) , k / 2) *
        (1 / alpha) * sqrt((1 + 1 / alpha) * Rcpp::trigamma( 1 + 1 / alpha) -1);

   return aux;
}


// [[Rcpp::export]]
NumericVector II_alpha (NumericVector alpha){
  NumericVector aux = (1 / alpha) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) - 1);
  return aux;
}


// [[Rcpp::export]]
NumericVector I_alpha (NumericVector alpha) {
  NumericVector aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) +
      pow(1 + Rcpp::digamma(1 + 1 / alpha), 2)  - 1);
  return aux;
}


// [[Rcpp::export]]
NumericVector prior_alpha(NumericVector alpha, int  k, int prior) {
  if (prior == 1){
    NumericVector aux = J_alpha(alpha, k);
    return aux;
  }
  if (prior == 2){
    NumericVector aux = II_alpha(alpha);
    return aux;
  }

  if (prior == 3){
    NumericVector aux = I_alpha(alpha);
    return aux;
  }
  return NULL;
}

// ##### PRIOR FOR (BETA,SIGMA2,ALPHA) (LOG-EXPONENTIAL POWER MODEL ONLY) ######

// [[Rcpp::export]]
NumericVector prior_LEP(NumericVector beta, float sigma2, NumericVector alpha,
                        int prior, bool logs ) {
  if (logs == false) {
    NumericVector aux = prior_LN(beta, sigma2, prior, logs = FALSE) *
      prior_alpha(alpha, beta.size(), prior);
    return(aux);
  } else{
    NumericVector aux = prior_LN(beta, sigma2, prior, logs = TRUE) +
      log(prior_alpha(alpha, beta.size(), prior));
    return(aux);
  }
}

/*
##### SAMPLING FROM A GENERALIZED INVERSE GAUSSIAN DISTRIBUTION
##### (REQUIRED FOR SEVERAL .LLOG FUNCTIONS)
##### BASED ON HOLMES AND HELD (2006)*/

// [[Rcpp::export]]
double r_GIG(double r) {
  double y = R::rnorm(0, 1);
  y = y * y;
  y = 1.0 + ((y - sqrt(y * (4.0 * r + y))) / (2.0 * r));
  double u = R::runif(0, 1);

    if (u <=  1 / (1 + y)) {
      double aux = r / y;
      return aux;
    } else {
      double aux = r * y;
      return aux;
    }
}


// [[Rcpp::export]]
double Log_aux(NumericVector lambda, double y, int j_nu, NumericVector nu,
        int prior){
double l_aux = sum(Rcpp::dgamma(lambda, y/2, 2/y, true)) -
    sum(Rcpp::dgamma(lambda, nu[j_nu]/2, 2/nu[j_nu], true)) +
      log(prior_nu(y, prior)[1] / prior_nu(nu[j_nu], prior)[1]);
  return l_aux;
}

// DENSITY FUNCTION OF A TRUNCATED EXPONENTIAL DISTRIBUTION
// (REQUIRED FOR BF.u.i.LEP ONLY)
// [[Rcpp::export]]
double d_texp(double x, double trunc) {
  if (x >= trunc) {
    double aux = exp(- (x - trunc));
    return(aux);
  } else {
    double aux = 0;
    return(aux);
  }
}

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
  }
  if (logs == true) {
    NumericVector aux = prior_LN(beta, sigma2, prior,
                                 true) + log(prior_nu(nu, prior));
    return aux;
  }

  return NULL;
}




// [[Rcpp::export]]
NumericVector J_alpha(NumericVector alpha, int k){
  NumericVector aux = pow(alpha * (alpha - 1) *

    Rcpp::trigamma(1.0 - 1.0 / alpha) /
      Rcpp::trigamma (1.0 / alpha) , k / 2) *
                                  (1 / alpha) *
                      sqrt((1 + 1 / alpha) * Rcpp::trigamma( 1 + 1 / alpha) -1);


   return aux;
}


// [[Rcpp::export]]
NumericVector II_alpha (NumericVector alpha){
  NumericVector aux = (1 / alpha) * sqrt((1 + 1 / alpha) * Rcpp::trigamma(1 + 1 / alpha) - 1);
  return aux;
}


// [[Rcpp::export]]
NumericVector I_alpha (NumericVector alpha) {
  NumericVector aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) + pow(1 + Rcpp::digamma(1 + 1 / alpha), 2)  - 1);
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
  }
  if (logs == true) {
    NumericVector aux = prior_LN(beta, sigma2, prior, logs = TRUE) +
      log(prior_alpha(alpha, beta.size(), prior));
    return(aux);
  }
  return(NULL);
}


// [[Rcpp::export]]
NumericVector r_GIG(double r,  int n = 1) {
  NumericVector y = Rcpp::rnorm(n);
  y = y * y;
  y = 1.0 + ((y - sqrt(y * (4.0 * r + y))) / (2.0 * r));
  NumericVector u = Rcpp::runif(n);
  NumericVector aux(n);
  for( int i = 0 ; i < n ; ++i) {
    if (u[i] <=  1 / (1 + y[i])) {
      aux[i] = r / y[i];
    } else {
      aux[i] = r * y[i];
    }
  }
  return aux;
}




// [[Rcpp::export]]
double Log_aux(NumericVector lambda, double y, int j_nu, NumericVector nu,
        int prior){
double l_aux = sum(Rcpp::dgamma(lambda, y/2, 2/y, true)) -
    sum(Rcpp::dgamma(lambda, nu[j_nu]/2, 2/nu[j_nu], true)) +
                                                   log(prior_nu(y, prior)[1] /
                                                     prior_nu(nu[j_nu], prior)[1]);
  return l_aux;
}






//
// double log_lik_LLOG ( NumericVector Time, NumericVector Cens, arma::mat X,
//                      arma::vec beta, double sigma2, int set, double eps_l,
//                      double eps_r) {
//   int n = Time.size();
//   NumericVector aux(n);
//
//   arma::vec MEAN = X * beta;
//   NumericVector MEANVec = NumericVector(MEAN.begin(), MEAN.end());
//
//
//   NumericVector sigma2Vec(n);
//   std::fill(sigma2Vec.begin(), sigma2Vec.end(), sigma2);
//
//   if (set == 1) {
//
//     NumericVector TimeLower(n);
//
//     for( int i = 0; i < n; i = i + 1){
//       if(Time[i] > eps_l){
//         TimeLower[i] = 1;
//       }  else{
//         TimeLower[i] = 0;
//       }
//     }
//
//     NumericVector aux = Cens * TimeLower *
//       log(Rcpp::plogis(log(Time + eps_r), MEANVec,
//                         sqrt(sigma2Vec), true, false) -
//                           R::plogis(log(abs(Time - eps_l)),
//                                        MEANVec,
//                                        sqrt(sigma2Vec))) +
//                                           (1 - TimeLower) *
//                                           log(R::plogis(log(Time + eps_r),
//                                                             MEAN,
//                                                             sqrt(sigma2Vec))) +
//                                                               (1 - Cens) *
//                                                               log(1 - Rcpp::plogis(log(Time),
//                                                                                     MEAN,
//                                                                                     sqrt(sigma2Vec)));
//   }
//   if (set == 0) {
//     NumericVector aux = Cens * (stats::dlogis(log(Time), MEAN, sqrt(sigma2),
//                                               true) -
//        log(Time)) + (1 - Cens) * log(1 - stats::plogis(log(Time), MEAN,
//            sqrt(sigma2)));
//   }
//   return(sum(aux));
// }



// double MH_nu_LST (int N, double omega2, NumericVector beta, double lambda, double nu0,
//                   int prior) {
//   int k = beta.size();
//   NumericVector ind(N);
//   NumericVector nu(N + 1);
//   nu[0] = nu0;
//
//   for (int j_nu = 0; j_nu < N; j_nu = j_nu + 1) {
//     double y = R::rnorm(nu[j_nu], sqrt(omega2));
//     bool ind_aux = y >= 0;
//     y = abs(y);
//     double u_aux = R::runif(0, 1);
//
//
//
//     double log_aux = sum(Rcpp::dgamma(lambda, shape = y/2, rate = y/2,
//                                  log = TRUE)) -
//                                    sum(stats::dgamma(lambda,
//                                                      shape = nu[j_nu] / 2,
//                                                      rate = nu[j_nu] / 2,
//                                                      log = TRUE)) +
//                                                        log(prior_nu(y, prior) /
//                                                          prior_nu(nu[j_nu], prior))

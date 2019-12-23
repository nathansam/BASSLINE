#include <Rmath.h>
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;


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
double J_alpha(double alpha, int k){
  double aux = pow(alpha * (alpha - 1) *
                              ::tgamma(1.0 - 1.0 / alpha) /
                                ::tgamma (1.0 / alpha) , k / 2) *
                                  (1 / alpha) *
                      sqrt((1 + 1 / alpha) * R::trigamma( 1 + 1 / alpha) -1);

  return aux;
}

// [[Rcpp::export]]
NumericVector IIJ_nu(NumericVector nu){
  NumericVector aux = sqrt(nu / (nu + 3)) *
          sqrt (Rcpp::trigamma(nu / 2) - Rcpp::trigamma((nu + 1) / 2) -
                 (2 * (nu + 3)) / (nu * pow(nu + 1,2)));
  return aux;
}



// [[Rcpp::export]]
NumericVector r_GIG(double r,  int n = 1) {
  NumericVector y = Rcpp::rnorm(n);
  y = pow(y, 2);
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
  return(aux);
}


// [[Rcpp::export]]
double II_alpha (double alpha){
  double aux = (1 / alpha) * sqrt((1 + 1 / alpha) * R::trigamma(1 + 1 / alpha) - 1);
  return(aux);
}


// [[Rcpp::export]]
double I_alpha (double alpha) {
  double aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
    R::trigamma(1 + 1 / alpha) + pow(1 + R::digamma(1 + 1 / alpha), 2)  - 1);
  return(aux);
}


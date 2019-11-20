#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}


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
NumericVector J_alpha(NumericVector alpha, int k){
  NumericVector aux = pow(alpha * (alpha - 1) *
                              Rcpp::gamma (1.0 - 1.0 / alpha) /
                                Rcpp::gamma (1.0 / alpha) , 2) *
                                  (1 / alpha) *
                      sqrt((1 + 1 / alpha) * Rcpp::trigamma( 1 + 1 / alpha) -1);

  return aux;
}

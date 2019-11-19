#include <Rcpp.h>

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
    float p = 1 + k / 2;
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

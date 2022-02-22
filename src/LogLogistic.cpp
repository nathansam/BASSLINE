#include "LogNormal.h"

// ######### SAMPLING FROM A GENERALIZED INVERSE GAUSSIAN DISTRIBUTION #########
//  BASED ON HOLMES AND HELD (2006)
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

// ##################### REJECTION SAMPLING FOR LAMBDA[i] ######################
// BASED ON HOLMES AND HELD (2006)
// [[Rcpp::export]]
Rcpp:: List RS_lambda_obs_LLOG(arma::vec logt, arma::mat X, arma::vec beta,
                               double sigma2, int obs, int N_AKS) {
  double lambda = 0;
  bool OK = false;
  int step = 0;

  double Z;
  double W;
  int n_AKS;

  while (OK == false) {
    step = step + 1;


    lambda = r_GIG((abs(logt[obs - 1] - X.row(obs - 1) * beta)[0]) /
      sqrt(sigma2));

    if (lambda != 0 && lambda != INFINITY) {
      double U = R::runif(0, 1);
      if (lambda > 4 / 3) {
        Z = 1;
        W = exp(-0.5 * lambda);
        n_AKS = 0;
        while (n_AKS <= N_AKS) {
          n_AKS = n_AKS + 1;
          Z = Z - (pow(n_AKS + 1, 2)) * pow(W, pow(n_AKS + 1, 2) - 1);
          if (Z > U) {
            OK = true;
          }
          n_AKS = n_AKS + 1;
          Z = Z + pow(n_AKS + 1, 2) * pow(W, pow(n_AKS + 1, 2) - 1);
          if (Z < U) {
            OK = false;
          }
        }
      } else {
        double H = 0.5 * log(2) + 2.5 * log(M_PI) - 2.5 *
          log(lambda) - ((M_PI * M_PI) / (2 * lambda)) + 0.5 *
          lambda;
        double lU = log(U);
        Z = 1;
        W = exp(- (M_PI * M_PI) / (2 * lambda));
        double K = lambda / (M_PI * M_PI);
        n_AKS = 0;
        while (n_AKS <= N_AKS) {
          n_AKS = n_AKS + 1;
          Z = Z - K * pow(W, (pow(n_AKS, 2) - 1));
          double aux = H + log(Z);

          if ( aux != R_NaReal && aux != -INFINITY && aux == INFINITY) {
            OK = true;
          }
          if ( aux != R_NaReal && aux < INFINITY && aux > -INFINITY && aux > lU) {
            OK = true;
          }
          n_AKS = n_AKS + 1;
          Z = Z + (pow(n_AKS + 1, 2)) * pow(W, pow(n_AKS + 1, 2) - 1);
          aux = H + log(Z);
          if ( aux != R_NaReal && aux == -INFINITY) {
            OK = false;
          }
          if (aux != R_NaReal && aux > -INFINITY && aux < INFINITY && aux < lU) {
            OK = false;
          }
        }
      }
    }
  }

  return List::create(Rcpp::Named("lambda", lambda),
                      Rcpp::Named("steps", step));
}



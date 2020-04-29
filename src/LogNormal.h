#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using arma::vec;
using Rcpp::List;
using arma::mat;



// [[Rcpp::export]]
double prior_LN(NumericVector beta, double sigma2, int prior, bool logs);

// [[Rcpp::export]]
NumericVector J_alpha(NumericVector alpha, int k);

// [[Rcpp::export]]
double J_alpha_single(double alpha, int k);

// [[Rcpp::export]]
NumericVector II_alpha (NumericVector alpha);

// [[Rcpp::export]]
double II_alpha_single (double alpha);


// [[Rcpp::export]]
NumericVector I_alpha (NumericVector alpha);

// [[Rcpp::export]]
double I_alpha_single (double alpha);


// [[Rcpp::export]]
NumericVector prior_alpha(NumericVector alpha, int  k, int prior);

// [[Rcpp::export]]
double prior_alpha_single(double alpha, int  k, int prior);




// [[Rcpp::export]]
arma::colvec rtnorm(int n, arma::vec lower, arma::vec upper,
                    arma::vec mu, arma::vec sd);
// [[Rcpp::export]]
NumericVector mvrnormArma(int n, arma::vec mu, arma::mat Sigma);

// [[Rcpp::export]]
arma::vec mvrnormArma2(int n, arma::vec mu, arma::mat Sigma);

// [[Rcpp::export]]
arma::vec logt_update_SMLN (arma::vec Time, arma::vec Cens,
                            arma::mat X, arma::vec beta, double sigma2,
                            bool set, double eps_l, double eps_r);

// [[Rcpp::export]]
arma::mat MCMC_LN_CPP (int N, int thin, int burn, arma::vec Time,
                       arma::vec Cens, arma::mat X, arma::vec beta0,
                       double sigma20, int prior, bool set, double eps_l,
                       double eps_r);

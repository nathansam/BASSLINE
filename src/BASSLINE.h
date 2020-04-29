/* C++ implementation of BASSLINE
  
* Author: Nathan Constantine-Cooke
  
*/


#define BASSLINE_H                                                              
                                                               

#include <RcppArmadillo.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using arma::vec;
using Rcpp::List;
using arma::mat;

// [[Rcpp::export]]
double prior_LN(NumericVector beta, double sigma2, int prior, bool logs);

// [[Rcpp::export]]
double prior_nu_single(double nu, int prior);

// [[Rcpp::export]]
double prior_LST(NumericVector beta, double sigma2, double nu, int prior,
                 bool logs);

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
NumericVector prior_LEP(NumericVector beta, float sigma2, NumericVector alpha,
                        int prior, bool logs);
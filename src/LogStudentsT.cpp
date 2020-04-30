#include "LogNormal.h"

/*##############################################################################
#################################### Priors ####################################
##############################################################################*/

// ################################ Prior for nu ###############################
// [[Rcpp::export]]
double prior_nu_single(double nu, int prior){
  double aux;

  if (prior == 2){
    aux = sqrt(nu / (nu + 3)) *
      sqrt (R::trigamma(nu / 2) - R::trigamma((nu + 1) / 2) -
      (2 * (nu + 3)) / (nu * pow(nu + 1, 2)));
    return aux;
  } else{
    return 0;
  }
}

// ######################## PRIOR FOR (BETA,SIGMA2,NU) #########################
// [[Rcpp::export]]
double prior_LST(NumericVector beta, double sigma2, double nu,
                 int prior, bool logs) {
  double aux;

  if (logs == false) {
    aux = prior_LN(beta, sigma2, prior, false) * prior_nu_single(nu, prior);
  } else {
    aux = prior_LN(beta, sigma2, prior, true) + log(prior_nu_single(nu, prior));
  }
  return aux;
}

// [[Rcpp::export]]
double alpha_nu(double nu0, double nu1, NumericVector lambda, int prior) {
  if (nu1 <= 0) {
    double aux = 0.0;
    return aux;
  } else {

    double aux = exp(sum(Rcpp::dgamma(lambda, nu1 / 2,
                                      2 / nu1, true))
                       - sum(Rcpp::dgamma(lambda,
                                          nu0 / 2,
                                          2 / nu0,
                                          true))) *
                                            (prior_nu_single(nu1, prior) /
                                                   prior_nu_single(nu0, prior));
    if (aux > 1){
      aux = 1.0;
    }
    return aux;
  }
}

// ###################### METROPOLIS-HASTING UPDATE OF NU ######################
// [[Rcpp::export]]
Rcpp::List MH_nu_LST(double omega2, NumericVector beta,
                     NumericVector lambda, double nu0, int prior) {


  int ind, ind_aux, aux;
  double nu, u_aux, log_aux;

  double y = R::rnorm(nu0, sqrt(omega2));

  if (y >= 0){
    ind_aux = 1;
  } else{
    ind_aux = 0;
  }

  y = abs(y);
  u_aux = R::runif(0, 1);

  log_aux = sum(Rcpp::dgamma(lambda, y/2, 2/y, true)) -
    sum(Rcpp::dgamma(lambda, nu0/2, 2/nu0, true))  +
    log(prior_nu_single(y, prior) /
      prior_nu_single(nu0, prior));


  if (log(u_aux) < log_aux){
    aux = 1;
  } else{
    aux = 0;
  }

  aux = aux * ind_aux;
  nu = aux * y + (1 - aux) * nu0;
  ind = aux;

  return List::create(Rcpp::Named("nu", nu), Rcpp::Named("ind", ind));
}


// LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL LST FUNCTIONS)
// [[Rcpp::export]]
double log_lik_LST(NumericVector Time, NumericVector Cens, arma::mat X,
                   arma::vec beta, double sigma2, double nu, bool set,
                   double eps_l, double eps_r) {
  const unsigned int n = Time.length();
  NumericVector aux(n);
  NumericVector MEAN = Rcpp::wrap(X * beta);

  NumericVector sigma2Vec(n, sigma2);
  NumericVector nuVec(n, nu);

  NumericVector TimeGreater(n);

  for (unsigned int i = 0; i < n; ++i){
    if (Time[i] > eps_l){
      TimeGreater[i] = 1;
    } else{
      TimeGreater[i] = 0;
    }
  }


  if (set == true) {
    aux = Cens * (TimeGreater *
      log(Rcpp::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2), nu) -
      Rcpp::pt((log(abs(Time - eps_l)) - MEAN) / sqrt(sigma2),
               nu)) +
                 (1 - TimeGreater) *
                 log(Rcpp::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2),
                              nu) - 0)) +
                                (1 - Cens) *
                                log(1 - Rcpp::pt((log(Time) - MEAN) /
                                                             sqrt(sigma2),
                                                 nu));
  } else {
    // set == 0
    aux = Cens * (Rcpp::dt((log(Time) - MEAN) / sqrt(sigma2),
                           nu, true) -
                             log(sqrt(sigma2) * Time)) +
                             (1 - Cens) * log(1 - Rcpp::pt((log(Time) - MEAN)
                                                             / sqrt(sigma2),
                                                             nu));
  }
  return sum(aux);
}


// ####################  MARGINAL POSTERIOR OF LAMBDA[obs] #####################
// [[Rcpp::export]]
double Post_lambda_obs_LST(const unsigned int obs, double  ref, arma::mat X,
                           arma::mat chain) {
  const unsigned int N = chain.n_rows;
  const unsigned int n = X.n_rows;
  const unsigned int k = X.n_cols;
  NumericVector aux1(N);
  NumericVector aux2(N);

  for (unsigned int iter = 0; iter < N; iter = iter + 1) {
    aux1[iter] = pow((chain(iter, obs + k + 1 + n) - X.row(obs - 1) *
      chain(iter, arma::span(0, k - 1)).t())[0], 2) / (chain(iter, k )) +
      chain(iter, k + 1);
    aux2[iter] = R::dgamma(ref, (chain(iter, (k + 1)) + 1) / 2,
                           2 / aux1[iter], false);
  }

  double aux = mean(aux2);
  return aux;
}

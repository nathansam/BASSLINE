#include "LogNormal.h"

/*##############################################################################
#################################### Priors ####################################
##############################################################################*/

// ###################### PRIOR FOR (BETA,SIGMA2,ALPHA)  ######################
// [[Rcpp::export]]
 NumericVector prior_LEP(NumericVector beta, float sigma2, NumericVector alpha,
                         int prior, bool logs ) {
   if (logs == false) {
     NumericVector aux = prior_LN(beta, sigma2, prior, logs = FALSE) *
       prior_alpha(alpha, beta.size(), prior);
     return aux;
   } else{
     NumericVector aux = prior_LN(beta, sigma2, prior, logs = TRUE) +
       log(prior_alpha(alpha, beta.size(), prior));
     return aux;
   }
 }


// ########## DENSITY FUNCTION OF A TRUNCATED EXPONENTIAL DISTRIBUTION #########
// [[Rcpp::export]]
double d_texp(double x, double trunc) {
  double aux;
  if (x >= trunc) {
    aux = exp(trunc - x);
    return aux;
  } else {
    aux = 0;
    return aux;
  }
}


// ################## METROPOLIS-HASTINMCMC UPDATE OF SIGMA2 ###################
// [[Rcpp::export]]
Rcpp::List MH_marginal_sigma2 (double omega2, arma::vec logt,
                               arma::mat X, arma::vec beta,
                               double alpha, double sigma20, int prior) {
  const unsigned int k = beta.n_elem;
  const unsigned int n = logt.n_elem;

  double p;
  int ind;

  if (prior == 1) {
    p = 1 + k / 2;
  } else {
    p = 1;
  }

  double sigma2;
  double  y_aux = R::rnorm(sigma20, sqrt(omega2));
  if (y_aux <= 0) {

    ind = 2;

  } else {
    double l1 = - ((n / 2) + p) * log(y_aux / sigma20);

    double l2 = - sum(pow(abs(logt - X * beta), alpha)) *
      (pow(y_aux , -alpha / 2) - pow(sigma20, -alpha / 2));

    double aux = l1 + l2;

    double u_aux = R::runif(0, 1);

    // THE FOLLOWING THREE LINES AVOID CRUSHES
    if (aux == R_NaReal) {
      sigma2 = sigma20;
      ind = 0;
      Rcpp::Rcout << "NA.sigma2\n";
    } else if (aux == INFINITY) {
      sigma2 = sigma20;
      ind = 0;
      Rcpp::Rcout << "-Inf.sigma2\n";
    } else if (aux == -INFINITY) {
      sigma2 = y_aux;
      ind = 1;
      Rcpp::Rcout << "Inf.sigma2\n";
    } else if (aux > -INFINITY && log(u_aux) < aux) {
      sigma2 = y_aux;
      ind = 1;
    } else{
      //aux > -INFINITY && log(u_aux) >= aux)
      sigma2 = sigma20;
      ind = 0;
    }
  }

  return List::create(Rcpp::Named("sigma2",sigma2), Rcpp::Named("ind",ind));
}

// ##### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF BETA[j] #####
// [[Rcpp::export]]
double alpha_beta(arma::vec beta_0, arma::vec beta_1,
                  arma::vec logt, arma::mat X, double sigma2, double alpha) {
  double l1 = -sum(pow(abs((logt - X * beta_1) / sqrt(sigma2)), alpha));
  double l2 = -sum(pow(abs((logt - X * beta_0) / sqrt(sigma2)), alpha));


  double aux = exp(l1 - l2);

  if (1 < aux){
    aux = 1;
  }

  return aux;
}


// ##### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF SIGMA2 ######
// [[Rcpp::export]]
double alpha_sigma2 (double sigma2_0, double sigma2_1, arma::vec logt,
                     arma::mat X, arma::vec beta, double alpha, int prior) {

  double aux;

  if (sigma2_1 <= 0) {
    aux = 0;
  } else {
    const unsigned int k = X.n_cols;
    const unsigned int n = X.n_rows;

    float p;
    if (prior == 1) {
      p = 1 + k / 2;
    } else {
      p = 1;
    }

    double l1 = - ((n / 2) + p) * log(sigma2_1 / sigma2_0);
    double l2 = - sum(pow(abs(logt - X * (beta)),alpha)) *
      (pow(sigma2_1, (-alpha / 2)) - pow(sigma2_0, (-alpha / 2)));

    aux = exp(l1 + l2);

    if (1 < aux ){
      aux = 1;
    }
  }

  return aux;
}


// ################### METROPOLIS-HASTINMCMC UPDATE OF ALPHA ###################
// [[Rcpp::export]]
Rcpp::List MH_marginal_alpha (double omega2, arma::vec logt,
                              arma::mat X, arma::vec beta, double sigma2,
                              double alpha0, int prior) {
  const unsigned int k = X.n_cols;
  const unsigned int n = X.n_rows;
  double alpha, y_aux, u_aux, aux, l0, l1;
  int ind;


  y_aux = R::rnorm(alpha0, omega2);
  if (y_aux >= 2 || y_aux <= 1) {
    alpha = alpha0;
    ind = 0;
  } else {

    l1 = n * log(y_aux / std::tgamma(1 / y_aux)) -
      sum(pow(abs((logt - X * beta) / sqrt(sigma2)) , y_aux));

    l0 = n * log(alpha0 / std::tgamma(1 / alpha0)) -
      sum(pow(abs((logt - X * beta) / sqrt(sigma2)), alpha0));

    aux = l1 - l0 + log(prior_alpha_single(y_aux, k, prior) /
      prior_alpha_single(alpha0, k, prior));

    u_aux = R::runif(0, 1);

    // THE FOLLOWING THREE LINES AVOID CRUSHES
    if (aux == R_NaN) {
      alpha = alpha0;
      ind = 0;
      Rcpp::Rcout << "NA.alpha\n";
    } else if (aux == -INFINITY) {
      alpha = alpha0;
      ind = 0;
      Rcpp::Rcout << "-Inf.alpha\n";
    } else if ( aux == INFINITY) {
      alpha = y_aux;
      ind = 1;
      Rcpp::Rcout << "Inf.alpha\n";
    } else if (log(u_aux) < aux) {
      alpha = y_aux;
      ind = 1;
    } else {
      alpha = alpha0;
      ind = 0;
    }
  }
  return List::create(Rcpp::Named("alpha", alpha), Rcpp::Named("ind",ind));
}




// ################## METROPOLIS-HASTINMCMC UPDATE OF BETA[j] ##################
// [[Rcpp::export]]
Rcpp::List MH_marginal_beta_j(double omega2, arma::vec logt,
                              arma::mat X, double sigma2, double alpha,
                              arma::vec beta0, int j) {

  int ind;
  vec beta;

  arma::vec y_aux = beta0;

  y_aux[j - 1] = R::rnorm(beta0(j - 1), sqrt(omega2));

  double aux = -sum(pow(abs((logt - (X * y_aux)) / sqrt(sigma2)), alpha)) +
    sum(pow(abs((logt - (X * beta0)) / sqrt(sigma2)) , alpha));
  double u_aux = R::runif(0, 1);


  if (aux == R_NaReal) {
    beta = beta0;
    ind = 0;
    Rcpp::Rcout << "NA.beta\n";
  } else if ( aux == -INFINITY) {
    beta = beta0;
    ind = 0;
    Rcpp::Rcout << "-Inf.beta\n";
  } else if (aux == INFINITY) {
    beta = y_aux;
    ind = 1;
    Rcpp::Rcout << "Inf.beta\n";
  } else if (log(u_aux) < aux) {
    beta = y_aux;
    ind = 1;
  } else{
    // (log(u_aux) >= aux)
    beta = beta0;
    ind = 0;
  }
  return List::create(Rcpp::Named("beta", beta), Rcpp::Named("ind",ind));
}


// ###### ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF ALPHA ######
// [[Rcpp::export]]
double alpha_alpha(double alpha0, double alpha1, arma::vec logt,
                   arma::mat X, arma::vec beta, double sigma2, int prior) {
  const unsigned int k = X.n_cols;
  const unsigned int n = X.n_rows;

  double aux;

  if (alpha1 < 1) {
    aux =  0;
  } else if (alpha1 > 2){
    aux = 0;
  } else {
    double l1 = n * log(alpha1 / std::tgamma(1 / alpha1)) -
      pow((1 / sqrt(sigma2)), alpha1) * sum(pow(abs(logt - X *
      beta), alpha1));
    double l0 = n * log(alpha0 / std::tgamma(1 / alpha0)) -
      (pow(1 / sqrt(sigma2), alpha0)) * sum(pow(abs(logt - X *
      beta) , alpha0));


    aux = exp(l1 - l0) * prior_alpha_single(alpha1, k, prior) /
      prior_alpha_single(alpha0, k, prior);

    if (aux > 1)
      aux = 1;
  }
  return aux;
}



// ########## DENSITY FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION ###########
// (BASED ON library normalp).
// BUGGED Currently broken on travis
// [[Rcpp::export]]
NumericVector d_normp(NumericVector x, NumericVector mu, NumericVector sigmap,
                      NumericVector p, bool logs = false) {

  const unsigned int n = x.length();

  NumericVector cost(n);
  NumericVector expon1(n);
  NumericVector expon2(n);
  NumericVector dsty(n);

  for (unsigned int i = 0; i < n; ++i){
    cost[i] = 2 * pow(p[i], 1 / p[i]) * std::tgamma(1 + 1 / p[i]) * sigmap[i];
    expon1[i] = pow(abs(x[i] - mu[i]), p[i]);
    expon2[i] = p[i] * pow(sigmap[i], p[i]);
    dsty[i] = (1 / cost[i]) * exp(-1 * expon1[i] / expon2[i]);
  }


  if (logs == true){
    dsty = log(dsty);
  }
  return dsty;
}


// ######## DISTRIBUTION FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION ########
// (BASED ON library normalp)
// NEEDS to be reworked to work for different argument types
// BUGGED Currently broken on travis
// [[Rcpp::export]]
NumericVector p_normp(NumericVector q, NumericVector mu,
                      NumericVector sigmap,
                      NumericVector p, bool lower_tail = true,
                      bool log_pr = false) {

  NumericVector z = (q - mu) / sigmap;
  const unsigned int n = p.length();
  NumericVector zz(n);
  NumericVector zp(n);
  NumericVector lzp(n);

  for( unsigned int i = 0; i < n; ++i){
    zz[i] = pow(abs(z[i]), p[i]);
    zp[i] = R::pgamma(zz[i], 1 / p[i], p[i], lower_tail, false);
    lzp[i] = R::pgamma(zz[i], 1 / p[i], p[i], lower_tail, true);

    if (z[i] < 0){
      zp[i] = 0.5 - exp(lzp[i] - log(2));
    } else {
      zp[i] = 0.5 + exp(lzp[i] - log(2));
    }

    if (log_pr == true){
      if (z[i] < 0){
        zp[i] = log(0.5 - exp(lzp[i] - log(2)));
      } else {
        zp[i] = log(0.5 + exp(lzp[i] - log(2)));
      }
    }
  }
  return zp;
}



// ########################## LOG-LIKELIHOOD FUNCTION ##########################
// [[Rcpp::export]]
double log_lik_LEP(NumericVector Time, NumericVector Cens, arma::mat X,
                   arma::vec beta, double sigma2, double alpha, bool set,
                   double eps_l, double eps_r) {
  const unsigned int n = Time.length();

  NumericVector aux(n);


  NumericVector MEAN = Rcpp::wrap(X * beta);

  NumericVector sigma2vec(n, sigma2);
  NumericVector alphavec(n, alpha);

  NumericVector SP(n);

  for (unsigned int i = 0; i < n; ++i){
    SP[i] = sqrt(sigma2vec[i]) *  pow(1 / alphavec[i], 1 / alphavec[i]);
  }

  if (set == true) {

    NumericVector TimeGreater(n);

    for (unsigned int i = 0; i < n; ++i){
      if (Time[i] > eps_l){
        TimeGreater[i] = 1;
      } else{
        TimeGreater[i] = 0;
      }
    }

    NumericVector aux1 = TimeGreater *
                     log(p_normp(log(Time + eps_r), MEAN, SP, alphavec) -
                       p_normp(log(abs(Time - eps_l)), MEAN, SP, alphavec)) +
                         (1 - TimeGreater) * log(p_normp(log(Time + eps_r),
                                                         MEAN, SP, alphavec));

    NumericVector aux2 = log(1 - p_normp(log(Time), MEAN, SP, alphavec));

    for (unsigned int i = 0; i < n; ++i){
      if (Cens[i] == 1){
        aux[i] = aux1[i];
      } else{
        aux[i] = aux2[i];
      }
    }

  }
  if (set == false) {
    NumericVector aux = Cens * (d_normp(log(Time), MEAN, SP,
                                        alphavec) - log(Time)) + (1 - Cens) *
                                          log(1 - p_normp(log(Time), MEAN, SP,
                                                          alphavec));
  }
  return sum(aux);
}

// ######################## MARGINAL POSTERIOR OF u[obs] #######################
// BUGGED
// [[Rcpp::export]]
double Post_u_obs_LEP(const unsigned int obs, double ref, arma::mat X,
                      arma::mat chain) {
  const unsigned int N = chain.n_rows;
  const unsigned int n = X.n_rows;
  const unsigned int k = X.n_cols;
  NumericVector aux1(N);

  for (unsigned int iter = 0; iter < N;  ++iter) {

    double trunc_aux = pow(abs((chain(iter, obs + k + 1 + n) - X.row(obs -1) *
                           chain(iter, arma::span(0, k - 1)).t())[0]) / sqrt(chain(iter, k)),
                           chain(iter, k + 1));

    aux1[iter] = d_texp(ref, trunc_aux);
  }
  double aux = mean(aux1);
  return aux;
}

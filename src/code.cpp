#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <RcppArmadillo.h>

using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::NumericMatrix;

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
  return 1;
}

// Non vectorised version of prior_nu
// Used for alpha_nu
// [[Rcpp::export]]
double prior_nu_single(double nu, int prior){

  if (prior == 2){
    double aux = sqrt(nu / (nu + 3)) *
      sqrt (R::trigamma(nu / 2) - R::trigamma((nu + 1) / 2) -
      (2 * (nu + 3)) / (nu * pow(nu + 1,2)));
    return aux;
  } else{
    return 0;
  }
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
double J_alpha_single(double alpha, int k){
  double aux = pow(alpha * (alpha - 1) *
    R::trigamma(1.0 - 1.0 / alpha) /
      R::trigamma (1.0 / alpha) , k / 2) *
        (1 / alpha) * sqrt((1 + 1 / alpha) * R::trigamma( 1 + 1 / alpha) -1);

  return aux;
}


// [[Rcpp::export]]
NumericVector II_alpha (NumericVector alpha){
  NumericVector aux = (1 / alpha) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) - 1);
  return aux;
}

// [[Rcpp::export]]
double II_alpha_single (double alpha){
  double aux = (1 / alpha) * sqrt((1 + 1 / alpha) *
    R::trigamma(1 + 1 / alpha) - 1);
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
double I_alpha_single (double alpha) {
  double aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
    R::trigamma(1 + 1 / alpha) +
    pow(1 + R::digamma(1 + 1 / alpha), 2)  - 1);
  return aux;
}


// [[Rcpp::export]]
NumericVector prior_alpha(NumericVector alpha, int  k, int prior) {
  NumericVector aux;
  if (prior == 1){
    aux = J_alpha(alpha, k);
    return aux;
  }
  if (prior == 2){
    aux = II_alpha(alpha);
    return aux;
  }

  if (prior == 3){
    aux = I_alpha(alpha);
    return aux;
  }
  return 1;
}

// [[Rcpp::export]]
double prior_alpha_single(double alpha, int  k, int prior) {

  double aux;

  if (prior == 1){
     aux = J_alpha_single(alpha, k);
    return aux;
  }
  if (prior == 2){
    aux = II_alpha_single(alpha);
    return aux;
  }

  if (prior == 3){
    aux = I_alpha_single(alpha);
    return aux;
  }
  return 0;
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


// DENSITY FUNCTION OF A TRUNCATED EXPONENTIAL DISTRIBUTION
// (REQUIRED FOR BF.u.i.LEP ONLY)
// [[Rcpp::export]]
double d_texp(double x, double trunc) {
  if (x >= trunc) {
    double aux = exp(- (x - trunc));
    return(aux);
  } else {
    double aux = 0;
    return aux;
  }
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
               (prior_nu_single(nu1, prior) / prior_nu_single(nu0, prior));
    if (aux > 1){
      aux = 1.0;
    }
  return aux;
  }
}


// METROPOLIS-HASTINMCMC UPDATE OF SIGMA2
// (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
// [[Rcpp::export]]
Rcpp::List MH_marginal_sigma2 (int N, double omega2, arma::vec logt,
                               arma::mat X, arma::vec beta,
                               double alpha, double sigma20, int prior) {
  const unsigned int k = beta.n_elem;
  const unsigned int n = logt.n_elem;

  double p;

  if (prior == 1) {
     p = 1 + k / 2;
  }
  if (prior == 2 || prior == 3) {
     p = 1;
  }

    NumericVector sigma2(N + 1);
    sigma2[0] = sigma20;
    NumericVector ind(N + 1);
    for (int i_s = 0; i_s < N; i_s = i_s + 1) {
      double  y_aux = R::rnorm(sigma2[i_s],sqrt(omega2));
      if (y_aux <= 0) {
        sigma2[i_s + 1] = sigma2[i_s];
        ind[i_s + 1] = 2;
      } else {
        double l1 = - ((n / 2) + p) * log(y_aux / sigma2[i_s]);

        double l2 = - sum(pow(abs(logt - X * beta), alpha)) *
          (pow(y_aux , -alpha / 2) - pow(sigma2[i_s], -alpha / 2));

        double aux = l1 + l2;

        double u_aux = R::runif(0, 1);

// THE FOLLOWING THREE LINES AVOID CRUSHES
        if (aux == R_NaN) {
          sigma2[i_s + 1] = sigma2[i_s];
          ind[i_s + 1] = 0;
          Rcpp::Rcout << "NA.sigma2\n";
        } else if (aux == INFINITY) {
          sigma2[i_s + 1] = sigma2[i_s];
          ind[i_s + 1] = 0;
          Rcpp::Rcout << "-Inf.sigma2\n";
        } else if (aux == -INFINITY) {
          sigma2[i_s + 1] = y_aux;
          ind[i_s + 1] = 1;
          Rcpp::Rcout << "Inf.sigma2\n";
        } else if (aux > -INFINITY && log(u_aux) < aux) {
          sigma2[i_s + 1] = y_aux;
          ind[i_s + 1] = 1;
        } else if (aux > -INFINITY && log(u_aux) >= aux) {
          sigma2[i_s + 1] = sigma2[i_s];
          ind[i_s + 1] = 0;
        }
      }
    }
    sigma2.erase(0);
    ind.erase(0);
    return List::create(Rcpp::Named("sigma2",sigma2), Rcpp::Named("ind",ind));
}

//METROPOLIS-HASTINMCMC UPDATE OF NU (REQUIRED FOR SEVERAL .LST FUNCTIONS)
// POSSIBLE
// [[Rcpp::export]]
Rcpp::List MH_nu_LST( unsigned int N, double omega2, NumericVector beta,
                      NumericVector lambda, double nu0, int prior) {
  NumericVector ind(N + 1);
  NumericVector nu(N + 1);
  nu[0] = nu0;

  for (int j_nu = 0; j_nu < N; j_nu = j_nu + 1) {
    double y = R::rnorm(nu[j_nu], sqrt(omega2));

    int ind_aux;
    if (y >= 0){
      ind_aux = 1;
    } else{
      ind_aux = 0;
    }

    y = abs(y);
    double u_aux = R::runif(0, 1);

    double log_aux = sum(Rcpp::dgamma(lambda, y/2, 2/y, true)) -
               sum(Rcpp::dgamma(lambda, nu[j_nu]/2, 2/nu[j_nu], true))  +
                   log(prior_nu_single(y, prior) /
                        prior_nu_single(nu[j_nu], prior));

    double aux;
    if (log(u_aux) < log_aux){
      aux = 1;
    } else{
      aux = 0;
    }

      aux = aux * ind_aux;
      nu[j_nu + 1] = aux * y + (1 - aux) * nu[j_nu];
      ind[j_nu + 1] = aux;
  }
  nu.erase(0);
  ind.erase(0);
  return List::create(Rcpp::Named("nu",nu), Rcpp::Named("ind",ind));
}

// ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF BETA[j]
//(REQUIRED FOR ML.LEP FUNCTION ONLY)
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


// ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF SIGMA2
// (REQUIRED FOR ML.LEP FUNCTION ONLY)
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
    }
    if (prior == 2) {
      p = 1;
    }
    if (prior == 3) {
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


// METROPOLIS-HASTINMCMC UPDATE OF ALPHA (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
// [[Rcpp::export]]
Rcpp::List MH_marginal_alpha (unsigned int N, double omega2, arma::vec logt,
                              arma::mat X, arma::vec beta, double sigma2,
                              double alpha0, int prior) {
  const unsigned int k = X.n_cols;
  const unsigned int n = X.n_rows;
  NumericVector alpha(N + 1);
  alpha[0] = alpha0;
  NumericVector ind(N + 1);

  for (int i_a = 0; i_a < N; i_a = i_a + 1) {
    double y_aux = R::rnorm(alpha[i_a], omega2);
    if (y_aux >= 2 || y_aux <= 1) {
      alpha[i_a + 1] = alpha[i_a];
      ind[i_a + 1] = 0;
    } else {

      double l1 = n * log(y_aux / std::tgamma(1 / y_aux)) -
        sum(pow(abs((logt - X * beta) / sqrt(sigma2)) , y_aux));

      double l0 = n * log(alpha[i_a] / std::tgamma(1 / alpha[i_a])) -
        sum(pow(abs((logt - X * beta) / sqrt(sigma2)), alpha[i_a]));

      double aux = l1 - l0 + log(prior_alpha_single(y_aux, k, prior) /
        prior_alpha_single(alpha[i_a], k, prior));

      double u_aux = R::runif(0, 1);

// THE FOLLOWING THREE LINES AVOID CRUSHES
      if (aux == R_NaN) {
        alpha[i_a + 1] = alpha[i_a];
        ind[i_a + 1] = 0;
        Rcpp::Rcout << "NA.alpha\n";
      } else if (aux == -INFINITY) {
        alpha[i_a + 1] = alpha[i_a];
        ind[i_a + 1] = 0;
        Rcpp::Rcout << "-Inf.alpha\n";
      } else if ( aux == INFINITY) {
        alpha[i_a + 1] = y_aux;
        ind[i_a + 1] = 1;
        Rcpp::Rcout << "Inf.alpha\n";
      } else if (log(u_aux) < aux) {
        alpha[i_a + 1] = y_aux;
        ind[i_a + 1] = 1;
      } else {
        alpha[i_a + 1] = alpha[i_a];
        ind[i_a + 1] = 0;
      }
    }
  }

  alpha.erase(0);
  ind.erase(0);

  return List::create(Rcpp::Named("alpha", alpha), Rcpp::Named("ind",ind));
}



// DISTRIBUTION FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
// (BASED ON library normalp)
// NEED TO FIX FOR WHEN p IS NumericVector
// [[Rcpp::export]]
NumericVector pnormp_cpp(NumericVector q, NumericVector mu,
                         NumericVector sigmap,
                         double p, bool lower_tail = true,
                         bool log_pr = false) {

      NumericVector z = (q - mu) / sigmap;
      NumericVector zz = pow(abs(z), p);
      NumericVector zp = Rcpp::pgamma(zz, 1 / p, p);
      NumericVector lzp = Rcpp::pgamma(zz, 1 / p, p, true);
      zp = ifelse(z < 0, 0.5 - exp(lzp - log(2)), 0.5 + exp(lzp - log(2)));
      if (log_pr == TRUE)
        zp = ifelse(z < 0, log(0.5 - exp(lzp - log(2))),
                     log(0.5 + exp(lzp - log(2))));
        return zp;
}


// METROPOLIS-HASTINMCMC UPDATE OF BETA[j]
// (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
// NEED TO FIX!!!!!!!!!
// [[Rcpp::export]]
Rcpp::List MH_marginal_beta_j(unsigned int N, double omega2, arma::vec logt,
                              arma::mat X, double sigma2, double alpha,
                              arma::vec beta0, int j) {
  const unsigned int k = beta0.n_elem;
  arma::mat beta = arma::zeros(k, N + 1);

  beta.col(0) = beta0;

  NumericVector ind(N + 1);

  for (int i_b = 0; i_b < N; i_b = i_b + 1){

    arma::vec y_aux = beta.col(i_b);

    y_aux[j - 1] = R::rnorm(beta(j, i_b), sqrt(omega2));

    double aux = -sum(pow(abs((logt - (X * y_aux)) / sqrt(sigma2)), alpha)) +
      sum(pow(abs((logt - (X * beta.col(i_b))) / sqrt(sigma2)) , alpha));
    double u_aux = R::runif(0, 1);


    if (aux == R_NaN) {
      beta.col(i_b + 1) = beta.col(i_b);
      ind[i_b + 1] = 0;
      Rcpp::Rcout << "NA.beta\n";
    } else if ( aux == -INFINITY) {
      beta.col(i_b + 1) = beta.col(i_b);
      ind[i_b + 1] = 0;
      Rcpp::Rcout << "-Inf.beta\n";
    } else if (aux == INFINITY) {
      beta.col(i_b + 1) = y_aux;
      ind[i_b + 1] = 1;
      Rcpp::Rcout << "Inf.beta\n";
    } else if (log(u_aux) < aux) {
      beta.col(i_b + 1) = y_aux;
      ind[i_b + 1] = 1;
    } else if (log(u_aux) >= aux) {
      beta.col(i_b + 1) = beta.col(i_b);
      ind[i_b + 1] = 0;
    }
  }
  //beta.shed_col(0);
  ind.erase(0);

  return List::create(Rcpp::Named("beta", beta.t()), Rcpp::Named("ind",ind));
}


// ACEPTANCE PROBABILITY FOR METROPOLIS-HASTINMCMC UPDATE OF ALPHA
// (REQUIRED FOR ML.LEP FUNCTION ONLY)
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




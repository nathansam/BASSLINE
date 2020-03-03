#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using arma::vec;
using Rcpp::List;
using arma::mat;

/*##############################################################################
#####################################PRIORS#####################################
##############################################################################*/

//########################### PRIOR FOR (BETA,SIGMA2) ##########################

// [[Rcpp::export]]
double prior_LN(NumericVector beta, double sigma2, int prior, bool logs){

  float p;
  double aux;

  int k = beta.length();

  if (prior == 1){
    p = 1.0 + k / 2.0;
  } else{
    p = 1;
  }

  if (logs == false){
    aux = pow(sigma2, -1 * p);
  }

  if (logs == true){
    aux = log(sigma2) * -p;
  }
  return aux;
}

// Prior for nu
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


//########### PRIOR FOR (BETA,SIGMA2,NU) (LOG-STUDENT'S T MODEL ONLY)###########

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



// BUGGED
// [[Rcpp::export]]
NumericVector J_alpha(NumericVector alpha, int k){
  NumericVector aux = pow(alpha * (alpha - 1.0) *
    Rcpp::gamma(1.0 - 1.0 / alpha) /
      Rcpp::gamma (1.0 / alpha) , k / 2.0) *
        (1.0 / alpha) *
          sqrt((1.0 + 1.0 / alpha) * Rcpp::trigamma( 1.0 + 1.0 / alpha) - 1.0);

   return aux;
}

// BUGGED
// [[Rcpp::export]]
double J_alpha_single(double alpha, int k){
  double aux = pow(alpha * (alpha - 1.0) * std::tgamma(1.0 - 1.0 / alpha) /
                   std::tgamma(1.0 / alpha) ,
                   (k / 2.0)) * (1.0 / alpha) *
                     sqrt((1.0 + 1.0 / alpha) *
                        R::trigamma(1.0 + 1.0 / alpha) - 1.0);
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

  double aux;
  if (x >= trunc) {
    aux = exp(trunc - x);
    return aux;
  } else {
    aux = 0;
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

//METROPOLIS-HASTINMCMC UPDATE OF NU (REQUIRED FOR SEVERAL .LST FUNCTIONS)
// POSSIBLE
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


// METROPOLIS-HASTINMCMC UPDATE OF ALPHA (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
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




// METROPOLIS-HASTINMCMC UPDATE OF BETA[j]
// (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
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

// DENSITY FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
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


// DISTRIBUTION FUNCTION OF THE EXPONENTIAL POWER DISTRIBUTION
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



// LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL .LEP FUNCTIONS)
// [[Rcpp::export]]
double log_lik_LEP(NumericVector Time, NumericVector Cens, arma::mat X,
                   arma::vec beta, double sigma2, double alpha, int set,
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

  if (set == 1) {

    NumericVector TimeGreater(n);

    for (unsigned int i = 0; i < n; ++i){
      if (Time[i] > eps_l){
        TimeGreater[i] = 1;
      } else{
        TimeGreater[i] = 0;
      }
    }

    NumericVector aux1 = TimeGreater * log(p_normp(log(Time + eps_r), MEAN, SP,
                                                   alphavec) -
                            p_normp(log(abs(Time - eps_l)),
                                    MEAN, SP, alphavec)) +
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
  if (set == 0) {
    NumericVector aux = Cens * (d_normp(log(Time), MEAN, SP,
                           alphavec) - log(Time)) + (1 - Cens) *
                             log(1 - p_normp(log(Time), MEAN, SP, alphavec));
  }
  return sum(aux);
}


// REJECTION SAMPLING FOR LAMBDA[i]
// (REQUIRED FOR SEVERAL .LLOG FUNCTIONS)
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
        double H = 0.5 * log(2) + 2.5 * log(PI) - 2.5 *
          log(lambda) - ((PI * PI) / (2 * lambda)) + 0.5 *
          lambda;
        double lU = log(U);
        Z = 1;
        W = exp(- (PI * PI) / (2 * lambda));
        double K = lambda / (PI * PI);
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

// MARGINAL POSTERIOR OF u[obs]
// (REQUIRED FOR BF.u.obs.LEP ONLY)
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

// MARGINAL POSTERIOR OF LAMBDA[obs] (REQUIRED FOR BF.lambda.obs.LST ONLY)
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


double rtnormsingle(double mu, double sd, double lower, double upper){
  if (lower == R_NegInf) {
    lower = -1e35;
  }
  if (upper == R_PosInf) {
    upper = 1e35;
  }


  double z, pz, u, slower, supper, tr, alpha;
  bool sample = 1;

  if(lower >= upper){
    return((lower + upper) / 2);
  }
  if(lower < -1e32 || upper > 1e32){

    if(lower < -1e32 && upper > 1e32){
      z = R::rnorm(mu, sd);
      return z;

    }else if(upper > 1e32){
        tr = (lower - mu) / sd;
      }else{
        tr = (mu - upper) / sd;
      }
      if(tr < 0){
        /* if sampling >0.5 of a normal density possibly quicker just to
           sample and reject */
        while(sample == 1){
          z = R::rnorm(0.0, 1.0);

          if(z > tr){
            sample = 0;
          }
        }
      }else{
        alpha = (tr + sqrt((tr * tr) + 4.0)) / 2.0;
        while(sample == 1){

          z = R::rexp(1.0/alpha) + tr;
          pz = - ((alpha - z) * (alpha - z) / 2.0);
          u = -R::rexp(1.0);
          if(u <= pz){
            sample = 0;
          }
        }
  }
  }else{

    slower = (lower - mu) / sd;
    supper = (upper - mu) / sd;

    while(sample == 1){
      z = R::runif(slower, supper);

      if(slower <= 0.0 && 0.0 <= supper){
        pz = -z * z / 2.0;
      }else{
        if(supper < 0.0){
          pz = (supper * supper - z * z) / 2.0;
        }else{
          pz = (slower * slower - z * z) / 2.0;
        }
      }

      u = - R::rexp(1.0);

      if(u < pz){
        sample = 0;
      }
    }
  }
  if(lower < -1e32){
    return(mu - z * sd);
  }else{
    return(z * sd + mu);
  }
}

// C++ adaptation of Jarrod Hadfield's MCMCglmm::rtnorm
// [[Rcpp::export]]
NumericVector rtnorm(int n, NumericVector lower, NumericVector upper,
                     NumericVector mu, NumericVector sd){

  if (mu.length() == 1) {
    mu = NumericVector (n, mu[0]);
  }

  if (sd.length() == 1) {
    sd = NumericVector (n, sd[0]);
  }

  if (lower.length() == 1) {
    lower = NumericVector (n, lower[0]);
  }

  if (upper.length() == 1) {
    upper = NumericVector (n, upper[0]);
  }

  NumericVector rv (n);
  for( int i = 0; i < n; i++){
    rv[i] = rtnormsingle(mu[i], sd[i], lower[i], upper[i]);
  }
  return rv;
}

// LOG-LIKELIHOOD FUNCTION (REQUIRED FOR SEVERAL LST FUNCTIONS)
// [[Rcpp::export]]
double log_lik_LST(NumericVector Time, NumericVector Cens, arma::mat X,
                   arma::vec beta, double sigma2, double nu, int set,
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


  if (set == 1) {
    aux = Cens * (TimeGreater *
      log(Rcpp::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2), nu) -
        Rcpp::pt((log(abs(Time - eps_l)) - MEAN) / sqrt(sigma2),
                 nu)) +
                   (1 - TimeGreater) *
                     log(Rcpp::pt((log(Time + eps_r) - MEAN) / sqrt(sigma2),
                                  nu) - 0)) +
                       (1 - Cens) *
                         log(1 - Rcpp::pt((log(Time) - MEAN) / sqrt(sigma2),
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
  return(sum(aux));
}

// [[Rcpp::export]]
NumericVector mvrnormArma(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  Y = arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
  Rcpp::NumericVector tmp = Rcpp::wrap(Y);
  tmp.attr("dim") = R_NilValue;

  return tmp;
}

// [[Rcpp::export]]
arma::vec mvrnormArma2(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  Y = arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
  return Y.t();
}

// [[Rcpp::export]]
NumericVector logt_update_SMLN (NumericVector Time, NumericVector Cens,
                                arma::mat X, arma::vec beta, double sigma2,
                                int set, double eps_l, double eps_r) {
  const int n = Time.length();
  NumericVector aux(n);

  NumericVector MEAN = Rcpp::wrap(X * beta);
  NumericVector sdVec (n, sqrt(sigma2));
  NumericVector maxUpper(n, 1e35);

  if (set == 1) {
    NumericVector TimeGreater(n);
    for (int i = 0; i < n; i++){
      if (Time[i] > eps_l){
        TimeGreater[i] = 1;
      }
    }

    NumericVector minLower(n, -1e35);

    aux = Cens * (TimeGreater *
      rtnorm(n, log(abs(Time - eps_l)), log(Time + eps_r), MEAN, sdVec) +
         (1 - TimeGreater) * rtnorm(n, minLower, log(Time + eps_r), MEAN,
           sdVec)) + (1 - Cens) * rtnorm(n, log(Time), maxUpper, MEAN, sdVec);
  } else{
    aux = Cens * log(Time) + (1 - Cens) *
      rtnorm(n, log(Time), maxUpper, MEAN, sdVec);
  }
  return(aux);
}



// [[Rcpp::export]]
arma::mat MCMC_LN_CPP (int N, int thin, int burn, arma::vec Time,
                           NumericVector Cens, arma::mat X, arma::vec beta0,
                           double sigma20, int prior, int set,
                           double eps_l, double eps_r) {

      int k = beta0.n_elem; // How many betas?
      int n = Time.n_elem; // How many observations?
      int N_aux = N / thin;

      double p;

      if (prior == 1) {
        p = 1 + k / 2;
      }
      if (prior == 2) {
        p = 1;
      }

      // empty matrix with N_aux + 1 rows and k cols
      arma::mat beta = arma::zeros(N_aux + 1, k);

      beta.row(0) = beta0.t();

      arma::vec sigma2(N_aux + 1);
      sigma2[0] = sigma20;

      arma::mat logt = arma::zeros(N_aux + 1, n);
      logt.row(0) = log(Time).t();

      arma::vec beta_aux = beta0;

      double sigma2_aux = sigma2[0];
      arma::vec logt_aux = log(Time);


      // Identity matrix
      arma::mat D = arma::eye<arma::mat>(X.n_cols, X.n_cols);

  for (int iter = 1; iter < N + 1; iter ++) {



    arma::vec  mu_aux = arma::solve(X.t() * X, D) * X.t() * logt_aux;
    arma::mat Sigma_aux = sigma2_aux * solve(X.t() * X, D);


    arma::vec beta_aux = mvrnormArma2(1, mu_aux, Sigma_aux);



    double shape_aux = (n + 2.0 * p - 2.0) / 2.0;

    arma::vec rate_aux = 0.5 * (logt_aux - X * beta_aux).t() * (logt_aux - X *
      beta_aux);
    
    
    int max_i = rate_aux.n_elem;

    for (int i = 0; i < max_i; i++){
      if (rate_aux[i] > 0 && std::isnan(rate_aux[i]) == false) {
        sigma2_aux = pow(R::rgamma(shape_aux, 1.0 / rate_aux[i]), -1);
      }

      

    }

    logt_aux = logt_update_SMLN(Rcpp::wrap(Time), Rcpp::wrap(Cens), X, beta_aux,
                                sigma2_aux, set, eps_l, eps_r);

      if (iter % thin == 0) {
        beta.row(iter / thin) = beta_aux.t();
        sigma2[iter / thin] = sigma2_aux;
        logt.row(iter / thin) = logt_aux.t();
      }

      if ((iter - 1) % 100000 == 0) {
        Rcpp::Rcout << "Iteration :" << iter << std::endl;
      }
  }



  arma::mat chain = arma::join_horiz(beta, sigma2);
  chain = arma::join_horiz(chain, logt);
  return chain;
}

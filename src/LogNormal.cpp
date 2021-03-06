#include "uses.h"

//########################### PRIOR FOR (BETA,SIGMA2) ##########################

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


// BUGGED
NumericVector J_alpha(NumericVector alpha, int k){
  NumericVector aux = pow(alpha * (alpha - 1.0) *
    Rcpp::gamma(1.0 - 1.0 / alpha) /
      Rcpp::gamma (1.0 / alpha) , k / 2.0) *
        (1.0 / alpha) *
        sqrt((1.0 + 1.0 / alpha) * Rcpp::trigamma( 1.0 + 1.0 / alpha) - 1.0);

  return aux;
}


// BUGGED
double J_alpha_single(double alpha, int k){
  double aux = pow(alpha * (alpha - 1.0) * std::tgamma(1.0 - 1.0 / alpha) /
                   std::tgamma(1.0 / alpha), (k / 2.0)) *
                     (1.0 / alpha) * sqrt((1.0 + 1.0 / alpha) *
                     R::trigamma(1.0 + 1.0 / alpha) - 1.0);
  return aux;
}


NumericVector II_alpha (NumericVector alpha){
  NumericVector aux = (1 / alpha) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) - 1);
  return aux;
}


double II_alpha_single (double alpha){
  double aux = (1 / alpha) * sqrt((1 + 1 / alpha) *
                R::trigamma(1 + 1 / alpha) - 1);
  return aux;
}


NumericVector I_alpha (NumericVector alpha) {
  NumericVector aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
    Rcpp::trigamma(1 + 1 / alpha) +
    pow(1 + Rcpp::digamma(1 + 1 / alpha), 2)  - 1);
  return aux;
}


double I_alpha_single (double alpha) {
  double aux = sqrt(1 / pow(alpha, 3)) * sqrt((1 + 1 / alpha) *
                    R::trigamma(1 + 1 / alpha) +
                    pow(1 + R::digamma(1 + 1 / alpha), 2)  - 1);
  return aux;
}


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
  return 0;
}


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

int checkInterrupt(int count){
  count ++;
  if (count == 100){
    Rcpp::checkUserInterrupt();
    count = 0;
  }
  return (count);
}

double rtnormsingle(double mu, double sd, double lower, double upper){
  // Which algorithm to use
  int alg;
  int count = 0;
  bool done = false;
  double y, a, u, rho;

  if (upper == R_PosInf){
    upper = 1e35;
  }

  if (lower == R_NegInf){
    lower = -1e35;
  }


  if (lower  > upper){
    // Bigger lower than upper is not valid
    alg = -1;

  } else if (((lower < 0) & (upper > 1e10)) |
    ((lower == -1e10) & (upper > 0)) |
    ((lower < 0) & (upper > 0) & (upper - lower > sqrt(2 * M_PI)))){
    /* standard "simulate from normal and reject if outside limits" method.
    Use if bounds are wide. */
    alg = 0;

  } else if (((lower >= 0) & (upper > lower + 2 * sqrt(exp(1)) /
    (lower + sqrt(lower * lower + 4)) *
      exp((lower * 2 - lower * sqrt(lower * lower + 4)) / 4)))){
    // rejection sampling with exponential proposal. Use if lower >> mean
    alg = 1;

  } else if ((upper <= 0) & (-lower > -upper + 2*sqrt(exp(1)) /
    (-upper + sqrt(upper * upper + 4)) * exp((upper * upper -
      -upper * sqrt(upper * upper + 4)) / 4))){
      // rejection sampling with exponential proposal. Use if upper << mean.
      alg = 2;
  }  else{
    /* rejection sampling with uniform proposal.
     * Use if bounds are narrow and central. */
    alg = 3;
  }

  if (alg == - 1){
    Rcpp::Rcout << "Lower limit of a truncation cannot be biggger than the "<<
      "upper limit" << std::endl;
    y = 0;
    return(0);
  }

  if (alg == 0){
    while (done == false){
      y = ::Rf_rnorm(0, 1);
      count = checkInterrupt(count);
      if ((y > lower) & (y < upper)){
        done = true;
      }
    }
  }

  if (alg == 1){
    while (done == false){
      a = (lower + sqrt(lower * lower + 4)) / 2;
      y = ::Rf_rexp(a) + lower;
      u = ::Rf_runif(0, 1);

      count = checkInterrupt(count);

      if ((u <= exp(-pow(y - a, 2) / 2)) & (y <= upper)){
        done = true;
      }
    }
  }

  if (alg == 2){
    while (done == false){
      a = (- upper + sqrt(upper * upper + 4)) / 2;
      y = ::Rf_rexp(a) - upper;
      u = ::Rf_runif(0, 1);

      count = checkInterrupt(count);

      if ((u <= exp(- pow ( y - a, 2) / 2)) & (y <= -lower)){
        done = true;
      }
    }
    y = - y;
  }

  if (alg == 3){
    while (done == false){
      y = ::Rf_runif(lower, upper);

      if (lower > 0){
        rho = exp((lower * lower - y * y) / 2);
      } else if (upper < 0){
        rho = exp ((upper * upper - y * y ) / 2);
      } else{
        rho = exp ((-y * y)/2);
      }

      u = Rf_runif(0, 1);

      count = checkInterrupt(count);

      if (u <= rho){
        done = true;
      }
    }

  }

  return (y * sd + mu);

}

arma::vec Vect(int n, arma::vec x){
  // If arma::vec is length 1 then repeat value in vector n times.
  if (x.n_elem == 1){
    double xVal = x[0];
    x = arma::vec(n);
    x.fill(xVal);
  }
  return x;
}

/* C++ version of truncated normal based on R function
 *  msm::rtnorm by Christopher Jackson */
arma::colvec rtnorm(int n, arma::vec lower, arma::vec upper,
                    arma::vec mu, arma::vec sd){

  mu = Vect(n, mu); sd = Vect(n, sd);
  lower = Vect(n, lower); upper = Vect(n, upper);

  // Algorithm works on a mean 0, sd 1 scale
  lower = (lower - mu) / sd;
  upper = (upper - mu) / sd;

  arma::colvec rv (n);
  for( int i = 0; i < n; i++){
    rv[i] = rtnormsingle(mu[i], sd[i], lower[i], upper[i]);
  }
  return rv;
}


NumericVector mvrnormArma(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  Y = arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
  Rcpp::NumericVector tmp = Rcpp::wrap(Y);
  tmp.attr("dim") = R_NilValue;

  return tmp;
}


arma::vec mvrnormArma2(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  Y = arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
  return Y.t();
}


arma::vec logt_update_SMLN (arma::vec Time, arma::vec Cens,
                            arma::mat X, arma::vec beta, double sigma2,
                            bool set, double eps_l, double eps_r) {
  int n = Time.n_elem;
  arma::vec aux(n);

  arma::vec MEAN = X * beta;
  arma::vec sdVec (n);
  sdVec.fill(sqrt(sigma2));

  arma::vec maxUpper (n);
  maxUpper.fill(1e35);
  arma::vec minLower(n);
  minLower.fill(-1e35);

  if (set == true) {
    arma::vec TimeGreater(n);

    for (int i = 0; i < n; i++){
      if (Time[i] > eps_l){
        TimeGreater[i] = 1;
      }
    }

    aux = Cens % (TimeGreater %
      rtnorm(n, log(abs(Time - eps_l)), log(Time + eps_r), MEAN, sdVec) +
      (1 - TimeGreater) % rtnorm(n, minLower, log(Time + eps_r), MEAN,
       sdVec)) + (1 - Cens) %
         rtnorm(n, log(Time), maxUpper, MEAN, sdVec);

  } else{
    aux = Cens % log(Time) + (1 - Cens) %
      rtnorm(n, log(Time), maxUpper, MEAN, sdVec);
  }


  return aux;
}


arma::mat MCMC_LN_CPP (int N, int thin, int burn, arma::vec Time,
                       arma::vec Cens, arma::mat X, arma::vec beta0,
                       double sigma20, int prior, bool set, double eps_l,
                       double eps_r) {

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

    logt_aux = logt_update_SMLN(Time, Cens, X, beta_aux,
                                sigma2_aux, set, eps_l, eps_r);

    if (iter % thin == 0) {
      beta.row(iter / thin) = beta_aux.t();
      sigma2[iter / thin] = sigma2_aux;
      logt.row(iter / thin) = logt_aux.t();
    }

    if ((iter) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    if ((iter) % 100000 == 0) {
      Rcpp::Rcout << "Iteration :" << iter << std::endl;
    }
  }
  arma::mat chain = arma::join_horiz(beta, sigma2);
  chain = arma::join_horiz(chain, logt);
  return chain;
}






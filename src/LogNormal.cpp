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


double rtnormsingle(double mu, double sd, double lower, double upper){

  double z, pz, u, slower, supper, tr, alpha;
  int count = 0; // Used to check for user interrupts in while loops
  bool sample = 1;

  if (lower == R_NegInf) {
    lower = -1e10;
  }
  if (upper == R_PosInf) {
    upper = 1e10;
  }

  if(lower >= upper){
    return((lower + upper) / 2);
  }
  if(lower < -1e9 || upper > 1e9){

    if(lower < -1e9 && upper > 1e9){
      z = ::Rf_rnorm(mu, sd);
      return z;

    }else if(upper > 1e9){
      tr = (lower - mu) / sd;
    }else{
      tr = (mu - upper) / sd;
    }
    if(tr < 0){
      /* if sampling >0.5 of a normal density possibly quicker just to
       sample and reject */

        while(sample == 1){

          count = count + 1;
          if (count == 100){
            Rcpp::checkUserInterrupt();
            count = 0;
          }

          z = ::Rf_rnorm(0.0, 1.0);

          if(z > tr){
            sample = 0;
          }
        }
    }else{
      alpha = (tr + sqrt((tr * tr) + 4.0)) / 2.0;
      while(sample == 1){
        count = count + 1;
        if (count == 100){
          Rcpp::checkUserInterrupt();
          count = 0;
        }

        z = ::Rf_rexp(1.0/alpha) + tr;
        pz = - ((alpha - z) * (alpha - z) / 2.0);
        u = -1 * ::Rf_rexp(1.0);
        if(u <= pz){
          sample = 0;
        }
      }
    }
  }else{

    slower = (lower - mu) / sd;
    supper = (upper - mu) / sd;

    while(sample == 1){
      count = count + 1;
      if (count == 100){
        Rcpp::checkUserInterrupt();
        count = 0;
      }

      z = ::Rf_runif(slower, supper);

      if(slower <= 0.0 && 0.0 <= supper){
        pz = -z * z / 2.0;
      }else{
        if(supper < 0.0){
          pz = (supper * supper - z * z) / 2.0;
        }else{
          pz = (slower * slower - z * z) / 2.0;
        }
      }

      u = -1 *  ::Rf_rexp(1.0);

      if(u < pz){
        sample = 0;
      }
    }
  }
  if(lower < -1e9){
    return(mu - z * sd);
  }else{
    return(z * sd + mu);
  }
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

// C++ adaptation of Jarrod Hadfield's MCMCglmm::rtnorm
arma::colvec rtnorm(int n, arma::vec lower, arma::vec upper,
                    arma::vec mu, arma::vec sd){

  mu = Vect(n, mu); sd = Vect(n, sd);
  lower = Vect(n, lower); upper = Vect(n, upper);

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
  maxUpper.fill(1e10);

  if (set == true) {
    arma::vec TimeGreater(n);

    for (int i = 0; i < n; i++){
      if (Time[i] > eps_l){
        TimeGreater[i] = 1;
      }
    }

    arma::vec minLower(n);
    minLower.fill(-1e10);

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
    // ISSUE HERE!
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






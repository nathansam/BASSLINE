// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// prior_LN
double prior_LN(NumericVector beta, double sigma2, int prior, bool logs);
RcppExport SEXP _BASSLINE_prior_LN(SEXP betaSEXP, SEXP sigma2SEXP, SEXP priorSEXP, SEXP logsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< bool >::type logs(logsSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_LN(beta, sigma2, prior, logs));
    return rcpp_result_gen;
END_RCPP
}
// prior_nu
NumericVector prior_nu(NumericVector nu, int prior);
RcppExport SEXP _BASSLINE_prior_nu(SEXP nuSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_nu(nu, prior));
    return rcpp_result_gen;
END_RCPP
}
// prior_nu_single
double prior_nu_single(double nu, int prior);
RcppExport SEXP _BASSLINE_prior_nu_single(SEXP nuSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_nu_single(nu, prior));
    return rcpp_result_gen;
END_RCPP
}
// prior_LST
NumericVector prior_LST(NumericVector beta, double sigma2, NumericVector nu, int prior, bool logs);
RcppExport SEXP _BASSLINE_prior_LST(SEXP betaSEXP, SEXP sigma2SEXP, SEXP nuSEXP, SEXP priorSEXP, SEXP logsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< bool >::type logs(logsSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_LST(beta, sigma2, nu, prior, logs));
    return rcpp_result_gen;
END_RCPP
}
// J_alpha
NumericVector J_alpha(NumericVector alpha, int k);
RcppExport SEXP _BASSLINE_J_alpha(SEXP alphaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(J_alpha(alpha, k));
    return rcpp_result_gen;
END_RCPP
}
// J_alpha_single
double J_alpha_single(double alpha, int k);
RcppExport SEXP _BASSLINE_J_alpha_single(SEXP alphaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(J_alpha_single(alpha, k));
    return rcpp_result_gen;
END_RCPP
}
// II_alpha
NumericVector II_alpha(NumericVector alpha);
RcppExport SEXP _BASSLINE_II_alpha(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(II_alpha(alpha));
    return rcpp_result_gen;
END_RCPP
}
// II_alpha_single
double II_alpha_single(double alpha);
RcppExport SEXP _BASSLINE_II_alpha_single(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(II_alpha_single(alpha));
    return rcpp_result_gen;
END_RCPP
}
// I_alpha
NumericVector I_alpha(NumericVector alpha);
RcppExport SEXP _BASSLINE_I_alpha(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(I_alpha(alpha));
    return rcpp_result_gen;
END_RCPP
}
// I_alpha_single
double I_alpha_single(double alpha);
RcppExport SEXP _BASSLINE_I_alpha_single(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(I_alpha_single(alpha));
    return rcpp_result_gen;
END_RCPP
}
// prior_alpha
NumericVector prior_alpha(NumericVector alpha, int k, int prior);
RcppExport SEXP _BASSLINE_prior_alpha(SEXP alphaSEXP, SEXP kSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_alpha(alpha, k, prior));
    return rcpp_result_gen;
END_RCPP
}
// prior_alpha_single
double prior_alpha_single(double alpha, int k, int prior);
RcppExport SEXP _BASSLINE_prior_alpha_single(SEXP alphaSEXP, SEXP kSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_alpha_single(alpha, k, prior));
    return rcpp_result_gen;
END_RCPP
}
// prior_LEP
NumericVector prior_LEP(NumericVector beta, float sigma2, NumericVector alpha, int prior, bool logs);
RcppExport SEXP _BASSLINE_prior_LEP(SEXP betaSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP priorSEXP, SEXP logsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< bool >::type logs(logsSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_LEP(beta, sigma2, alpha, prior, logs));
    return rcpp_result_gen;
END_RCPP
}
// r_GIG
double r_GIG(double r);
RcppExport SEXP _BASSLINE_r_GIG(SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(r_GIG(r));
    return rcpp_result_gen;
END_RCPP
}
// d_texp
double d_texp(double x, double trunc);
RcppExport SEXP _BASSLINE_d_texp(SEXP xSEXP, SEXP truncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type trunc(truncSEXP);
    rcpp_result_gen = Rcpp::wrap(d_texp(x, trunc));
    return rcpp_result_gen;
END_RCPP
}
// alpha_nu
double alpha_nu(double nu0, double nu1, NumericVector lambda, int prior);
RcppExport SEXP _BASSLINE_alpha_nu(SEXP nu0SEXP, SEXP nu1SEXP, SEXP lambdaSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< double >::type nu1(nu1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha_nu(nu0, nu1, lambda, prior));
    return rcpp_result_gen;
END_RCPP
}
// MH_marginal_sigma2
Rcpp::List MH_marginal_sigma2(int N, double omega2, arma::vec logt, arma::mat X, arma::vec beta, double alpha, double sigma20, int prior);
RcppExport SEXP _BASSLINE_MH_marginal_sigma2(SEXP NSEXP, SEXP omega2SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP sigma20SEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type omega2(omega2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma20(sigma20SEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_marginal_sigma2(N, omega2, logt, X, beta, alpha, sigma20, prior));
    return rcpp_result_gen;
END_RCPP
}
// MH_nu_LST
Rcpp::List MH_nu_LST(const unsigned int N, double omega2, NumericVector beta, NumericVector lambda, double nu0, int prior);
RcppExport SEXP _BASSLINE_MH_nu_LST(SEXP NSEXP, SEXP omega2SEXP, SEXP betaSEXP, SEXP lambdaSEXP, SEXP nu0SEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type omega2(omega2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_nu_LST(N, omega2, beta, lambda, nu0, prior));
    return rcpp_result_gen;
END_RCPP
}
// alpha_beta
double alpha_beta(arma::vec beta_0, arma::vec beta_1, arma::vec logt, arma::mat X, double sigma2, double alpha);
RcppExport SEXP _BASSLINE_alpha_beta(SEXP beta_0SEXP, SEXP beta_1SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP sigma2SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_1(beta_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha_beta(beta_0, beta_1, logt, X, sigma2, alpha));
    return rcpp_result_gen;
END_RCPP
}
// alpha_sigma2
double alpha_sigma2(double sigma2_0, double sigma2_1, arma::vec logt, arma::mat X, arma::vec beta, double alpha, int prior);
RcppExport SEXP _BASSLINE_alpha_sigma2(SEXP sigma2_0SEXP, SEXP sigma2_1SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma2_0(sigma2_0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_1(sigma2_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha_sigma2(sigma2_0, sigma2_1, logt, X, beta, alpha, prior));
    return rcpp_result_gen;
END_RCPP
}
// MH_marginal_alpha
Rcpp::List MH_marginal_alpha(unsigned int N, double omega2, arma::vec logt, arma::mat X, arma::vec beta, double sigma2, double alpha0, int prior);
RcppExport SEXP _BASSLINE_MH_marginal_alpha(SEXP NSEXP, SEXP omega2SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP alpha0SEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type omega2(omega2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_marginal_alpha(N, omega2, logt, X, beta, sigma2, alpha0, prior));
    return rcpp_result_gen;
END_RCPP
}
// MH_marginal_beta_j
Rcpp::List MH_marginal_beta_j(const unsigned int N, double omega2, arma::vec logt, arma::mat X, double sigma2, double alpha, arma::vec beta0, int j);
RcppExport SEXP _BASSLINE_MH_marginal_beta_j(SEXP NSEXP, SEXP omega2SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP beta0SEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type omega2(omega2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_marginal_beta_j(N, omega2, logt, X, sigma2, alpha, beta0, j));
    return rcpp_result_gen;
END_RCPP
}
// alpha_alpha
double alpha_alpha(double alpha0, double alpha1, arma::vec logt, arma::mat X, arma::vec beta, double sigma2, int prior);
RcppExport SEXP _BASSLINE_alpha_alpha(SEXP alpha0SEXP, SEXP alpha1SEXP, SEXP logtSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(alpha_alpha(alpha0, alpha1, logt, X, beta, sigma2, prior));
    return rcpp_result_gen;
END_RCPP
}
// d_normp
NumericVector d_normp(NumericVector x, NumericVector mu, NumericVector sigmap, NumericVector p, bool logs);
RcppExport SEXP _BASSLINE_d_normp(SEXP xSEXP, SEXP muSEXP, SEXP sigmapSEXP, SEXP pSEXP, SEXP logsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmap(sigmapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type logs(logsSEXP);
    rcpp_result_gen = Rcpp::wrap(d_normp(x, mu, sigmap, p, logs));
    return rcpp_result_gen;
END_RCPP
}
// p_normp
NumericVector p_normp(NumericVector q, NumericVector mu, NumericVector sigmap, NumericVector p, bool lower_tail, bool log_pr);
RcppExport SEXP _BASSLINE_p_normp(SEXP qSEXP, SEXP muSEXP, SEXP sigmapSEXP, SEXP pSEXP, SEXP lower_tailSEXP, SEXP log_prSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmap(sigmapSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< bool >::type log_pr(log_prSEXP);
    rcpp_result_gen = Rcpp::wrap(p_normp(q, mu, sigmap, p, lower_tail, log_pr));
    return rcpp_result_gen;
END_RCPP
}
// log_lik_LEP
double log_lik_LEP(NumericVector Time, NumericVector Cens, arma::mat X, arma::vec beta, double sigma2, double alpha, int set, double eps_l, double eps_r);
RcppExport SEXP _BASSLINE_log_lik_LEP(SEXP TimeSEXP, SEXP CensSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP setSEXP, SEXP eps_lSEXP, SEXP eps_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Cens(CensSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type set(setSEXP);
    Rcpp::traits::input_parameter< double >::type eps_l(eps_lSEXP);
    Rcpp::traits::input_parameter< double >::type eps_r(eps_rSEXP);
    rcpp_result_gen = Rcpp::wrap(log_lik_LEP(Time, Cens, X, beta, sigma2, alpha, set, eps_l, eps_r));
    return rcpp_result_gen;
END_RCPP
}
// RS_lambda_obs_LLOG
Rcpp:: List RS_lambda_obs_LLOG(arma::vec logt, arma::mat X, arma::vec beta, double sigma2, int obs, int N_AKS);
RcppExport SEXP _BASSLINE_RS_lambda_obs_LLOG(SEXP logtSEXP, SEXP XSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP obsSEXP, SEXP N_AKSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type N_AKS(N_AKSSEXP);
    rcpp_result_gen = Rcpp::wrap(RS_lambda_obs_LLOG(logt, X, beta, sigma2, obs, N_AKS));
    return rcpp_result_gen;
END_RCPP
}
// Post_u_obs_LEP
double Post_u_obs_LEP(const unsigned int obs, double ref, arma::mat X, arma::mat chain);
RcppExport SEXP _BASSLINE_Post_u_obs_LEP(SEXP obsSEXP, SEXP refSEXP, SEXP XSEXP, SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type ref(refSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type chain(chainSEXP);
    rcpp_result_gen = Rcpp::wrap(Post_u_obs_LEP(obs, ref, X, chain));
    return rcpp_result_gen;
END_RCPP
}
// Post_lambda_obs_LST
double Post_lambda_obs_LST(const unsigned int obs, double ref, arma::mat X, arma::mat chain);
RcppExport SEXP _BASSLINE_Post_lambda_obs_LST(SEXP obsSEXP, SEXP refSEXP, SEXP XSEXP, SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type ref(refSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type chain(chainSEXP);
    rcpp_result_gen = Rcpp::wrap(Post_lambda_obs_LST(obs, ref, X, chain));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BASSLINE_prior_LN", (DL_FUNC) &_BASSLINE_prior_LN, 4},
    {"_BASSLINE_prior_nu", (DL_FUNC) &_BASSLINE_prior_nu, 2},
    {"_BASSLINE_prior_nu_single", (DL_FUNC) &_BASSLINE_prior_nu_single, 2},
    {"_BASSLINE_prior_LST", (DL_FUNC) &_BASSLINE_prior_LST, 5},
    {"_BASSLINE_J_alpha", (DL_FUNC) &_BASSLINE_J_alpha, 2},
    {"_BASSLINE_J_alpha_single", (DL_FUNC) &_BASSLINE_J_alpha_single, 2},
    {"_BASSLINE_II_alpha", (DL_FUNC) &_BASSLINE_II_alpha, 1},
    {"_BASSLINE_II_alpha_single", (DL_FUNC) &_BASSLINE_II_alpha_single, 1},
    {"_BASSLINE_I_alpha", (DL_FUNC) &_BASSLINE_I_alpha, 1},
    {"_BASSLINE_I_alpha_single", (DL_FUNC) &_BASSLINE_I_alpha_single, 1},
    {"_BASSLINE_prior_alpha", (DL_FUNC) &_BASSLINE_prior_alpha, 3},
    {"_BASSLINE_prior_alpha_single", (DL_FUNC) &_BASSLINE_prior_alpha_single, 3},
    {"_BASSLINE_prior_LEP", (DL_FUNC) &_BASSLINE_prior_LEP, 5},
    {"_BASSLINE_r_GIG", (DL_FUNC) &_BASSLINE_r_GIG, 1},
    {"_BASSLINE_d_texp", (DL_FUNC) &_BASSLINE_d_texp, 2},
    {"_BASSLINE_alpha_nu", (DL_FUNC) &_BASSLINE_alpha_nu, 4},
    {"_BASSLINE_MH_marginal_sigma2", (DL_FUNC) &_BASSLINE_MH_marginal_sigma2, 8},
    {"_BASSLINE_MH_nu_LST", (DL_FUNC) &_BASSLINE_MH_nu_LST, 6},
    {"_BASSLINE_alpha_beta", (DL_FUNC) &_BASSLINE_alpha_beta, 6},
    {"_BASSLINE_alpha_sigma2", (DL_FUNC) &_BASSLINE_alpha_sigma2, 7},
    {"_BASSLINE_MH_marginal_alpha", (DL_FUNC) &_BASSLINE_MH_marginal_alpha, 8},
    {"_BASSLINE_MH_marginal_beta_j", (DL_FUNC) &_BASSLINE_MH_marginal_beta_j, 8},
    {"_BASSLINE_alpha_alpha", (DL_FUNC) &_BASSLINE_alpha_alpha, 7},
    {"_BASSLINE_d_normp", (DL_FUNC) &_BASSLINE_d_normp, 5},
    {"_BASSLINE_p_normp", (DL_FUNC) &_BASSLINE_p_normp, 6},
    {"_BASSLINE_log_lik_LEP", (DL_FUNC) &_BASSLINE_log_lik_LEP, 9},
    {"_BASSLINE_RS_lambda_obs_LLOG", (DL_FUNC) &_BASSLINE_RS_lambda_obs_LLOG, 6},
    {"_BASSLINE_Post_u_obs_LEP", (DL_FUNC) &_BASSLINE_Post_u_obs_LEP, 4},
    {"_BASSLINE_Post_lambda_obs_LST", (DL_FUNC) &_BASSLINE_Post_lambda_obs_LST, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_BASSLINE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

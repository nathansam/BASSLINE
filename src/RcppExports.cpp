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
// J_alpha
double J_alpha(double alpha, int k);
RcppExport SEXP _BASSLINE_J_alpha(SEXP alphaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(J_alpha(alpha, k));
    return rcpp_result_gen;
END_RCPP
}
// IIJ_nu
NumericVector IIJ_nu(NumericVector nu);
RcppExport SEXP _BASSLINE_IIJ_nu(SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(IIJ_nu(nu));
    return rcpp_result_gen;
END_RCPP
}
// r_GIG
NumericVector r_GIG(double r, int n);
RcppExport SEXP _BASSLINE_r_GIG(SEXP rSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(r_GIG(r, n));
    return rcpp_result_gen;
END_RCPP
}
// II_alpha
double II_alpha(double alpha);
RcppExport SEXP _BASSLINE_II_alpha(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(II_alpha(alpha));
    return rcpp_result_gen;
END_RCPP
}
// I_alpha
double I_alpha(double alpha);
RcppExport SEXP _BASSLINE_I_alpha(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(I_alpha(alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BASSLINE_prior_LN", (DL_FUNC) &_BASSLINE_prior_LN, 4},
    {"_BASSLINE_J_alpha", (DL_FUNC) &_BASSLINE_J_alpha, 2},
    {"_BASSLINE_IIJ_nu", (DL_FUNC) &_BASSLINE_IIJ_nu, 1},
    {"_BASSLINE_r_GIG", (DL_FUNC) &_BASSLINE_r_GIG, 2},
    {"_BASSLINE_II_alpha", (DL_FUNC) &_BASSLINE_II_alpha, 1},
    {"_BASSLINE_I_alpha", (DL_FUNC) &_BASSLINE_I_alpha, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BASSLINE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

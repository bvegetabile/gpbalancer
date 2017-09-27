// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// construct_sqexp
NumericMatrix construct_sqexp(NumericMatrix X, NumericMatrix M, double sig, double noise);
RcppExport SEXP _gpbalancer_construct_sqexp(SEXP XSEXP, SEXP MSEXP, SEXP sigSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_sqexp(X, M, sig, noise));
    return rcpp_result_gen;
END_RCPP
}
// construct_poly
NumericMatrix construct_poly(NumericMatrix X, double p, double sig0, double sig1, double noise);
RcppExport SEXP _gpbalancer_construct_poly(SEXP XSEXP, SEXP pSEXP, SEXP sig0SEXP, SEXP sig1SEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type sig0(sig0SEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_poly(X, p, sig0, sig1, noise));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _gpbalancer_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// par_sq_exp
arma::mat par_sq_exp(arma::mat design_x, arma::vec vec_theta, double sig_noise);
RcppExport SEXP _gpbalancer_par_sq_exp(SEXP design_xSEXP, SEXP vec_thetaSEXP, SEXP sig_noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type design_x(design_xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vec_theta(vec_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sig_noise(sig_noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(par_sq_exp(design_x, vec_theta, sig_noise));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _gpbalancer_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpbalancer_construct_sqexp", (DL_FUNC) &_gpbalancer_construct_sqexp, 4},
    {"_gpbalancer_construct_poly", (DL_FUNC) &_gpbalancer_construct_poly, 5},
    {"_gpbalancer_timesTwo", (DL_FUNC) &_gpbalancer_timesTwo, 1},
    {"_gpbalancer_par_sq_exp", (DL_FUNC) &_gpbalancer_par_sq_exp, 3},
    {"_gpbalancer_rcpp_hello", (DL_FUNC) &_gpbalancer_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpbalancer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

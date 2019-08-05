// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// run_hdbm_mcmc
arma::vec run_hdbm_mcmc(arma::vec Y, arma::vec A, arma::mat M, arma::mat C1, arma::mat C2, arma::vec beta_m, arma::vec alpha_a, double pi_m, double pi_a, int burnin, int nsamples);
RcppExport SEXP _hdbm_run_hdbm_mcmc(SEXP YSEXP, SEXP ASEXP, SEXP MSEXP, SEXP C1SEXP, SEXP C2SEXP, SEXP beta_mSEXP, SEXP alpha_aSEXP, SEXP pi_mSEXP, SEXP pi_aSEXP, SEXP burninSEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C1(C1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C2(C2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_m(beta_mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< double >::type pi_m(pi_mSEXP);
    Rcpp::traits::input_parameter< double >::type pi_a(pi_aSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(run_hdbm_mcmc(Y, A, M, C1, C2, beta_m, alpha_a, pi_m, pi_a, burnin, nsamples));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hdbm_run_hdbm_mcmc", (DL_FUNC) &_hdbm_run_hdbm_mcmc, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_hdbm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
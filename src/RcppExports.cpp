// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kosmic_impl
List kosmic_impl(NumericVector input_vector, int decimals, int bootstrap, double t1min, double t1max, double t2min, double t2max, double sd, double tol);
RcppExport SEXP _tidykosmic_kosmic_impl(SEXP input_vectorSEXP, SEXP decimalsSEXP, SEXP bootstrapSEXP, SEXP t1minSEXP, SEXP t1maxSEXP, SEXP t2minSEXP, SEXP t2maxSEXP, SEXP sdSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_vector(input_vectorSEXP);
    Rcpp::traits::input_parameter< int >::type decimals(decimalsSEXP);
    Rcpp::traits::input_parameter< int >::type bootstrap(bootstrapSEXP);
    Rcpp::traits::input_parameter< double >::type t1min(t1minSEXP);
    Rcpp::traits::input_parameter< double >::type t1max(t1maxSEXP);
    Rcpp::traits::input_parameter< double >::type t2min(t2minSEXP);
    Rcpp::traits::input_parameter< double >::type t2max(t2maxSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(kosmic_impl(input_vector, decimals, bootstrap, t1min, t1max, t2min, t2max, sd, tol));
    return rcpp_result_gen;
END_RCPP
}
// kosmic_resamples_impl
List kosmic_resamples_impl(NumericVector results, NumericVector counts, int replicates, NumericVector settings);
RcppExport SEXP _tidykosmic_kosmic_resamples_impl(SEXP resultsSEXP, SEXP countsSEXP, SEXP replicatesSEXP, SEXP settingsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type results(resultsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< int >::type replicates(replicatesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type settings(settingsSEXP);
    rcpp_result_gen = Rcpp::wrap(kosmic_resamples_impl(results, counts, replicates, settings));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tidykosmic_kosmic_impl", (DL_FUNC) &_tidykosmic_kosmic_impl, 9},
    {"_tidykosmic_kosmic_resamples_impl", (DL_FUNC) &_tidykosmic_kosmic_resamples_impl, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_tidykosmic(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

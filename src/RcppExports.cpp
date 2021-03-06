// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_dist_euclidean
NumericMatrix cpp_dist_euclidean(const NumericMatrix& inMat);
RcppExport SEXP _macam_cpp_dist_euclidean(SEXP inMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type inMat(inMatSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dist_euclidean(inMat));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dist_maxnorm
NumericMatrix cpp_dist_maxnorm(const NumericMatrix& inMat);
RcppExport SEXP _macam_cpp_dist_maxnorm(SEXP inMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type inMat(inMatSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dist_maxnorm(inMat));
    return rcpp_result_gen;
END_RCPP
}
// cpp_quantile
double cpp_quantile(const NumericVector& inVec, const double q);
RcppExport SEXP _macam_cpp_quantile(SEXP inVecSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type inVec(inVecSEXP);
    Rcpp::traits::input_parameter< const double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_quantile(inVec, q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_greater
IntegerMatrix cpp_greater(const NumericMatrix& inMat, const double inThreshold);
RcppExport SEXP _macam_cpp_greater(SEXP inMatSEXP, SEXP inThresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type inMat(inMatSEXP);
    Rcpp::traits::input_parameter< const double >::type inThreshold(inThresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_greater(inMat, inThreshold));
    return rcpp_result_gen;
END_RCPP
}
// cpp_distOneZero
IntegerMatrix cpp_distOneZero(const NumericMatrix& inMat, const std::string& inMethod, double s);
RcppExport SEXP _macam_cpp_distOneZero(SEXP inMatSEXP, SEXP inMethodSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type inMat(inMatSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type inMethod(inMethodSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_distOneZero(inMat, inMethod, s));
    return rcpp_result_gen;
END_RCPP
}
// cpp_twinSurrogate
NumericMatrix cpp_twinSurrogate(const NumericMatrix& original_e, int dim, int num_iter, double s, const std::string& surrogate_option, const std::string& initial_point, const std::string& distance_method, int point_per_year, const std::string& s_update, int n_twin_threshold, bool output_message);
RcppExport SEXP _macam_cpp_twinSurrogate(SEXP original_eSEXP, SEXP dimSEXP, SEXP num_iterSEXP, SEXP sSEXP, SEXP surrogate_optionSEXP, SEXP initial_pointSEXP, SEXP distance_methodSEXP, SEXP point_per_yearSEXP, SEXP s_updateSEXP, SEXP n_twin_thresholdSEXP, SEXP output_messageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type original_e(original_eSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type surrogate_option(surrogate_optionSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type initial_point(initial_pointSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type distance_method(distance_methodSEXP);
    Rcpp::traits::input_parameter< int >::type point_per_year(point_per_yearSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type s_update(s_updateSEXP);
    Rcpp::traits::input_parameter< int >::type n_twin_threshold(n_twin_thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type output_message(output_messageSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_twinSurrogate(original_e, dim, num_iter, s, surrogate_option, initial_point, distance_method, point_per_year, s_update, n_twin_threshold, output_message));
    return rcpp_result_gen;
END_RCPP
}
// cpp_nanMatrix
NumericMatrix cpp_nanMatrix(int nrow, int ncol);
RcppExport SEXP _macam_cpp_nanMatrix(SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_nanMatrix(nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_macam_cpp_dist_euclidean", (DL_FUNC) &_macam_cpp_dist_euclidean, 1},
    {"_macam_cpp_dist_maxnorm", (DL_FUNC) &_macam_cpp_dist_maxnorm, 1},
    {"_macam_cpp_quantile", (DL_FUNC) &_macam_cpp_quantile, 2},
    {"_macam_cpp_greater", (DL_FUNC) &_macam_cpp_greater, 2},
    {"_macam_cpp_distOneZero", (DL_FUNC) &_macam_cpp_distOneZero, 3},
    {"_macam_cpp_twinSurrogate", (DL_FUNC) &_macam_cpp_twinSurrogate, 11},
    {"_macam_cpp_nanMatrix", (DL_FUNC) &_macam_cpp_nanMatrix, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_macam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

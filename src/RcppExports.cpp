// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// RcppOwenT01
double RcppOwenT01(double h, double a, int jmax, double cutpoint);
RcppExport SEXP OwenQ_RcppOwenT01(SEXP hSEXP, SEXP aSEXP, SEXP jmaxSEXP, SEXP cutpointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type jmax(jmaxSEXP);
    Rcpp::traits::input_parameter< double >::type cutpoint(cutpointSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppOwenT01(h, a, jmax, cutpoint));
    return rcpp_result_gen;
END_RCPP
}
// RcppOwenT
double RcppOwenT(double h, double a, int jmax, double cutpoint);
RcppExport SEXP OwenQ_RcppOwenT(SEXP hSEXP, SEXP aSEXP, SEXP jmaxSEXP, SEXP cutpointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type jmax(jmaxSEXP);
    Rcpp::traits::input_parameter< double >::type cutpoint(cutpointSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppOwenT(h, a, jmax, cutpoint));
    return rcpp_result_gen;
END_RCPP
}
// RcppOwenQ1
NumericVector RcppOwenQ1(int nu, double t, NumericVector delta, NumericVector R, int jmax, double cutpoint);
RcppExport SEXP OwenQ_RcppOwenQ1(SEXP nuSEXP, SEXP tSEXP, SEXP deltaSEXP, SEXP RSEXP, SEXP jmaxSEXP, SEXP cutpointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type jmax(jmaxSEXP);
    Rcpp::traits::input_parameter< double >::type cutpoint(cutpointSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppOwenQ1(nu, t, delta, R, jmax, cutpoint));
    return rcpp_result_gen;
END_RCPP
}

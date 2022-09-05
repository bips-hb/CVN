// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fisherTestGreater
double fisherTestGreater(int a, int b, int c, int d);
RcppExport SEXP _CVN_fisherTestGreater(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(fisherTestGreater(a, b, c, d));
    return rcpp_result_gen;
END_RCPP
}
// midPFisherTestGreater
double midPFisherTestGreater(int a, int b, int c, int d);
RcppExport SEXP _CVN_midPFisherTestGreater(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(midPFisherTestGreater(a, b, c, d));
    return rcpp_result_gen;
END_RCPP
}
// convertRawReports2TablesRcpp
Rcpp::DataFrame convertRawReports2TablesRcpp(Rcpp::IntegerMatrix reports, int n_drugs, int n_events);
RcppExport SEXP _CVN_convertRawReports2TablesRcpp(SEXP reportsSEXP, SEXP n_drugsSEXP, SEXP n_eventsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type reports(reportsSEXP);
    Rcpp::traits::input_parameter< int >::type n_drugs(n_drugsSEXP);
    Rcpp::traits::input_parameter< int >::type n_events(n_eventsSEXP);
    rcpp_result_gen = Rcpp::wrap(convertRawReports2TablesRcpp(reports, n_drugs, n_events));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CVN_fisherTestGreater", (DL_FUNC) &_CVN_fisherTestGreater, 4},
    {"_CVN_midPFisherTestGreater", (DL_FUNC) &_CVN_midPFisherTestGreater, 4},
    {"_CVN_convertRawReports2TablesRcpp", (DL_FUNC) &_CVN_convertRawReports2TablesRcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CVN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

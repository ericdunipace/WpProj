// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "WpProj_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sufficientStatistics
Rcpp::List sufficientStatistics(const NumericMatrix& X_, const NumericMatrix& Y_, const NumericMatrix& theta_, const Rcpp::List& options_);
RcppExport SEXP _WpProj_sufficientStatistics(SEXP X_SEXP, SEXP Y_SEXP, SEXP theta_SEXP, SEXP options_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type options_(options_SEXP);
    rcpp_result_gen = Rcpp::wrap(sufficientStatistics(X_, Y_, theta_, options_));
    return rcpp_result_gen;
END_RCPP
}
// xtyUpdate
matrix xtyUpdate(const NumericMatrix& X_, const NumericMatrix& Y_, const NumericMatrix& theta_, const NumericVector& result_, const Rcpp::List& options_);
RcppExport SEXP _WpProj_xtyUpdate(SEXP X_SEXP, SEXP Y_SEXP, SEXP theta_SEXP, SEXP result_SEXP, SEXP options_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type result_(result_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type options_(options_SEXP);
    rcpp_result_gen = Rcpp::wrap(xtyUpdate(X_, Y_, theta_, result_, options_));
    return rcpp_result_gen;
END_RCPP
}
// W2penalized
SEXP W2penalized(SEXP X_, SEXP Y_, SEXP theta_, SEXP family_, SEXP penalty_, SEXP groups_, SEXP unique_groups_, SEXP group_weights_, SEXP lambda_, SEXP nlambda_, SEXP lmin_ratio_, SEXP alpha_, SEXP gamma_, SEXP tau_, SEXP scale_factor_, SEXP penalty_factor_, SEXP opts_);
RcppExport SEXP _WpProj_W2penalized(SEXP X_SEXP, SEXP Y_SEXP, SEXP theta_SEXP, SEXP family_SEXP, SEXP penalty_SEXP, SEXP groups_SEXP, SEXP unique_groups_SEXP, SEXP group_weights_SEXP, SEXP lambda_SEXP, SEXP nlambda_SEXP, SEXP lmin_ratio_SEXP, SEXP alpha_SEXP, SEXP gamma_SEXP, SEXP tau_SEXP, SEXP scale_factor_SEXP, SEXP penalty_factor_SEXP, SEXP opts_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type family_(family_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type penalty_(penalty_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type groups_(groups_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type unique_groups_(unique_groups_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type group_weights_(group_weights_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nlambda_(nlambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lmin_ratio_(lmin_ratio_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha_(alpha_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type tau_(tau_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type scale_factor_(scale_factor_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type penalty_factor_(penalty_factor_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts_(opts_SEXP);
    rcpp_result_gen = Rcpp::wrap(W2penalized(X_, Y_, theta_, family_, penalty_, groups_, unique_groups_, group_weights_, lambda_, nlambda_, lmin_ratio_, alpha_, gamma_, tau_, scale_factor_, penalty_factor_, opts_));
    return rcpp_result_gen;
END_RCPP
}
// cost_calculation_
Rcpp::NumericMatrix cost_calculation_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, const double p);
RcppExport SEXP _WpProj_cost_calculation_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_calculation_(A_, B_, p));
    return rcpp_result_gen;
END_RCPP
}
// pbClean
void pbClean();
RcppExport SEXP _WpProj_pbClean() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    pbClean();
    return R_NilValue;
END_RCPP
}
// test
void test();
RcppExport SEXP _WpProj_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test();
    return R_NilValue;
END_RCPP
}
// selVarMeanGen
Rcpp::NumericMatrix selVarMeanGen(const SEXP& X_, const SEXP& theta_, const SEXP& beta_);
RcppExport SEXP _WpProj_selVarMeanGen(SEXP X_SEXP, SEXP theta_SEXP, SEXP beta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type beta_(beta_SEXP);
    rcpp_result_gen = Rcpp::wrap(selVarMeanGen(X_, theta_, beta_));
    return rcpp_result_gen;
END_RCPP
}
// sinkhorn_
Rcpp::List sinkhorn_(Rcpp::NumericVector p_, Rcpp::NumericVector q_, Rcpp::NumericMatrix cost_matrix_, double epsilon, int niterations);
RcppExport SEXP _WpProj_sinkhorn_(SEXP p_SEXP, SEXP q_SEXP, SEXP cost_matrix_SEXP, SEXP epsilonSEXP, SEXP niterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q_(q_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type cost_matrix_(cost_matrix_SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type niterations(niterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sinkhorn_(p_, q_, cost_matrix_, epsilon, niterations));
    return rcpp_result_gen;
END_RCPP
}
// transport_C_
Rcpp::List transport_C_(const Rcpp::NumericVector& mass_a_, const Rcpp::NumericVector& mass_b_, const Rcpp::NumericMatrix& cost_matrix_, const Rcpp::CharacterVector& method_, double epsilon_, int niter_);
RcppExport SEXP _WpProj_transport_C_(SEXP mass_a_SEXP, SEXP mass_b_SEXP, SEXP cost_matrix_SEXP, SEXP method_SEXP, SEXP epsilon_SEXP, SEXP niter_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_a_(mass_a_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_b_(mass_b_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cost_matrix_(cost_matrix_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_(epsilon_SEXP);
    Rcpp::traits::input_parameter< int >::type niter_(niter_SEXP);
    rcpp_result_gen = Rcpp::wrap(transport_C_(mass_a_, mass_b_, cost_matrix_, method_, epsilon_, niter_));
    return rcpp_result_gen;
END_RCPP
}
// transport_
Rcpp::List transport_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, double p, double ground_p, const Rcpp::CharacterVector& method_, bool a_sort, double epsilon_, int niter_);
RcppExport SEXP _WpProj_transport_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP, SEXP ground_pSEXP, SEXP method_SEXP, SEXP a_sortSEXP, SEXP epsilon_SEXP, SEXP niter_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ground_p(ground_pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< bool >::type a_sort(a_sortSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_(epsilon_SEXP);
    Rcpp::traits::input_parameter< int >::type niter_(niter_SEXP);
    rcpp_result_gen = Rcpp::wrap(transport_(A_, B_, p, ground_p, method_, a_sort, epsilon_, niter_));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_
double wasserstein_(const Rcpp::NumericVector& mass_, const Rcpp::NumericMatrix& cost_, const double p, const Rcpp::IntegerVector& from_, const Rcpp::IntegerVector& to_);
RcppExport SEXP _WpProj_wasserstein_(SEXP mass_SEXP, SEXP cost_SEXP, SEXP pSEXP, SEXP from_SEXP, SEXP to_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_(mass_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cost_(cost_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type from_(from_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type to_(to_SEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_(mass_, cost_, p, from_, to_));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_p_iid_
double wasserstein_p_iid_(const SEXP& X_, const SEXP& Y_, double p);
RcppExport SEXP _WpProj_wasserstein_p_iid_(SEXP X_SEXP, SEXP Y_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_p_iid_(X_, Y_, p));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_p_iid_p_
double wasserstein_p_iid_p_(const SEXP& X_, const SEXP& Y_, double p);
RcppExport SEXP _WpProj_wasserstein_p_iid_p_(SEXP X_SEXP, SEXP Y_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_p_iid_p_(X_, Y_, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_WpProj_sufficientStatistics", (DL_FUNC) &_WpProj_sufficientStatistics, 4},
    {"_WpProj_xtyUpdate", (DL_FUNC) &_WpProj_xtyUpdate, 5},
    {"_WpProj_W2penalized", (DL_FUNC) &_WpProj_W2penalized, 17},
    {"_WpProj_cost_calculation_", (DL_FUNC) &_WpProj_cost_calculation_, 3},
    {"_WpProj_pbClean", (DL_FUNC) &_WpProj_pbClean, 0},
    {"_WpProj_test", (DL_FUNC) &_WpProj_test, 0},
    {"_WpProj_selVarMeanGen", (DL_FUNC) &_WpProj_selVarMeanGen, 3},
    {"_WpProj_sinkhorn_", (DL_FUNC) &_WpProj_sinkhorn_, 5},
    {"_WpProj_transport_C_", (DL_FUNC) &_WpProj_transport_C_, 6},
    {"_WpProj_transport_", (DL_FUNC) &_WpProj_transport_, 8},
    {"_WpProj_wasserstein_", (DL_FUNC) &_WpProj_wasserstein_, 5},
    {"_WpProj_wasserstein_p_iid_", (DL_FUNC) &_WpProj_wasserstein_p_iid_, 3},
    {"_WpProj_wasserstein_p_iid_p_", (DL_FUNC) &_WpProj_wasserstein_p_iid_p_, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_WpProj(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

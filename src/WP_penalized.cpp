#include "limbs_types.h"
#include <vector>
#include "utils.h"
#include "wp_solver.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include "eta_progress_bar.h"


using namespace Rcpp;


//[[Rcpp::export]]
SEXP WPpenalized(SEXP X_,
                 SEXP Y_,
                 SEXP theta_,
                 SEXP power_,
                 SEXP family_,
                 SEXP penalty_,
                 SEXP groups_,
                 SEXP unique_groups_,
                 SEXP group_weights_,
                 SEXP lambda_,
                 SEXP nlambda_,
                 SEXP lmin_ratio_,
                 SEXP alpha_,
                 SEXP gamma_,
                 SEXP tau_,
                 SEXP scale_factor_,
                 SEXP penalty_factor_,
                 SEXP opts_)
{
  
  
  const matMap X(as<matMap >(X_));
  // const matMap Y_copy(as<matMap >(Y_));
  // const matMap theta_copy(as<matMap >(theta_));
  
  matrix Y = Rcpp::as<matMap >(Y_);//Y_copy;
  matrix theta = Rcpp::as<matMap >(theta_); //theta_copy;
  
  double power = Rcpp::as<double>(power_);
  
  const int S = Y.cols();
  const int p = X.rows();
  const int N = X.cols();
  
  matrix mu =  theta.transpose() * X; //start with ''true'' mu
  
  matrix xtx = matrix::Zero(p, p);
  matrix xty = matrix::Zero(p, 1); //may be resized in suff stat function
  
  vector scale_factor(as<vector>(scale_factor_));
  const vectorI groups(as<vectorI>(groups_));
  const vectorI unique_groups(as<vectorI>(unique_groups_));
  
  
  // In glmnet, we minimize
  //   1/(2n) * ||y - X * beta||^2 + lambda * ||beta||_1
  // which is equivalent to minimizing
  //   1/2 * ||y - X * beta||^2 + n * lambda * ||beta||_1
  vector group_weights(as<vector>(group_weights_));
  
  
  std::vector<vector> lambda(as< std::vector<vector> >(lambda_));
  
  vector lambda_tmp;
  lambda_tmp = lambda[0];
  
  int nl = as<int>(nlambda_);
  vector lambda_base(nl);
  
  int nlambda = lambda_tmp.size();
  
  List opts(opts_);
  const int maxit        = as<int>(opts["maxit"]);
  const double tol       = as<double>(opts["tol"]);
  const double alpha     = as<double>(alpha_);
  const double gamma     = as<double>(gamma_);
  const double tau       = as<double>(tau_);
  const int infm_maxit   = as<int>(opts["infm_maxit"]);
  const bool display_progress = as<bool>(opts["display_progress"]);
  CharacterVector method(as<CharacterVector>(opts["method"]));
  std::string transport_method = as<std::string>(opts["transport_method"]);
  const int model_size   = as<int>(opts["model_size"]);
  const bool not_same    = as<bool>(opts["not_same"]);
  double epsilon    = as<double>(opts["epsilon"]);
  double OTmaxit    = as<int>(opts["OTmaxit"]);
  const bool same = !(not_same);
  bool selection = false;
  // const double pseudo_obs = as<double>(opts["pseudo_observations"]);
  
  CharacterVector family(as<CharacterVector>(family_));
  std::string penalty(as< std::string >(penalty_));
  vector penalty_factor(as<vector>(penalty_factor_));
  
  vector obs_weight(N * S);
  obs_weight.fill(1.0);
  
  matrix xty_old = xty;
  // matrix xtx_old = xtx;
  
  //order indices of x * theta
  matrixI idx_mu(S,N);
  
  
  //change scale factor to make estimation easier, let's say
  if ( scale_factor.size() == 0 ) {
    // if ( (scale_factor.size() == 0) && (penalty != "selection.lasso")  ) {
    scale_factor.resize(X.rows());
    scale_factor = X.rowwise().norm();
  } //else if ( (penalty == "selection.lasso")  && (scale_factor.size() != 0)) {
  
  // initialize pointers
  WpSolver *solver = NULL; // solver doesn't point to anything yet
  
  // WpSolver(  int s_, //diff
  //            const refMat  &X_,
  //            const refMat &Y_,
  //            const double &power_,
  //            const vectorI &groups_,
  //            const vectorI &unique_groups_,
  //            vector &group_weights_,
  //            vector &penalty_factor_,
  //            const vector &scale_factor_,
  //            const double tol_ = 1e-6)
  // initialize classes
  solver = new WpSolver(S, X, Y, power, groups, unique_groups,
                            group_weights, penalty_factor,
                            scale_factor, tol);
  
  // initialize wpsolve
  solver->init_wpsolve();
  
  // generate lambda vector
  double lmax = 0.0;
  lmax = solver->compute_lambda_zero(penalty); //
  
  bool provided_lambda = false;
  std::string elasticnettxt(".net");
  bool is_net_pen = penalty.find(elasticnettxt) != std::string::npos;
  
  if (nlambda < 1) {
    double lmin = as<double>(lmin_ratio_) * lmax;
    
    // lambda_base.setLinSpaced(nl, std::log(lmin), std::log(lmax));
    lambda_base.setLinSpaced(nl, std::log(lmax), std::log(lmin));
    lambda_base = lambda_base.array().exp();
    // lambda_base(0) = 0.0;
    // lambda_base(nl-1) = 0.0;
    nlambda = lambda_base.size();
    
    lambda_tmp.resize(nlambda);
    
    if (is_net_pen)
    {
      lambda_tmp = (lambda_base.array() / alpha).matrix(); // * n; //
    } else
    {
      lambda_tmp = lambda_base; // * n; //
    }
    
  } else {
    provided_lambda = true;
    lambda_tmp = lambda[0];
  }
  
  // results matrix
  matrix beta = matrix::Zero(p * S, nlambda);
  
  // diagnostic checks
  IntegerVector niter( nlambda );
  double ilambda = 0.0;
  // double num_tol = Eigen::NumTraits<double>::dummy_precision();
  
  
  // progress bar
  if(display_progress) {
    Rcpp::Rcout << "\n";
  }
  ETAProgressBar pb;
  Progress prog( nlambda, display_progress, pb );
  
  
  // if (method(0) == "projection") {
  //   lambda_tmp.reverseInPlace();
  // }
  
  // const matrix xty_original = xty;
  // vector best_lambda(p);
  
  
  for (int i = 0; i < nlambda; i++)
  {
    // vectors to save current and last value of coefficients

    if (i % 10 == 0) {
    Rcpp::checkUserInterrupt();
    }
    
    //set current lambda
    ilambda = lambda_tmp(i);
    // Rcpp::Rcout << ilambda <<"\n";
    if(i == 0) {
      //intitalize with lambda and other parameters
      solver->init(ilambda, penalty,
                   alpha, gamma, tau);
      // if(method(0) != "projection") solver->beta_ones();
    } else {
      solver->init_warm(ilambda);
    }
    
    niter(i)     = solver->solve(maxit);
    
    
    // save coefficient otherwise
    beta.col(i)  = solver->get_beta();
    
    // break if larger than max coef
    if( countNonZero(beta.col(i)) > model_size) {
      beta.conservativeResize(Eigen::NoChange, i+1);
      lambda_tmp.conservativeResize(i+1);
      break;
    }
    
    //update progress bar
    // if(display_progress){
    prog.increment();
    // }
    // if (innerIter[i] > 1 && same && method(0) != "projection") {
    //   xty = xty_old;
    //   solver->init_warm_xty();
    // } //really slows things down!
  } //end loop over lambda values
  
  double d = solver->get_d();
  
  delete solver;
  solver = NULL;
  
  // if (method(0) != "projection") {
  //   lambda_tmp.reverseInPlace();
  // }
  
  return Rcpp::List::create(Named("beta")       = Rcpp::wrap(beta),
                            Named("lambda")     = lambda_tmp,
                            Named("niter")      = niter,
                            Named("d")          = d);
  
}


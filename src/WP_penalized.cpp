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
                 SEXP power_,
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
  
  matrix Y = Rcpp::as<matMap >(Y_);//Y_copy;

  double power = Rcpp::as<double>(power_);
  
  const int S = Y.cols();
  const int p = X.rows();
  const int N = X.cols();
  
  vector scale_factor(as<vector>(scale_factor_));
  const vectorI groups(as<vectorI>(groups_));
  const vectorI unique_groups(as<vectorI>(unique_groups_));
  
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
  const bool display_progress = as<bool>(opts["display_progress"]);
  const int model_size   = as<int>(opts["model_size"]);
  // const double pseudo_obs = as<double>(opts["pseudo_observations"]);
  
  std::string penalty(as< std::string >(penalty_));
  vector penalty_factor(as<vector>(penalty_factor_));
  
  vector obs_weight(N * S);
  obs_weight.fill(1.0);
  
  //order indices of x * theta
  matrixI idx_mu(S,N);
  
  
  //change scale factor to make estimation easier, let's say
  if ( scale_factor.size() == 0 ) {
    scale_factor.resize(X.rows());
    scale_factor = X.rowwise().norm();
  } 
  
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

  
  // progress bar
  if(display_progress) {
    Rcpp::Rcout << "\n";
  }
  ETAProgressBar pb;
  Progress prog( nlambda, display_progress, pb );
  
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
    if(i == 0) {
      //intitalize with lambda and other parameters
      solver->init(ilambda, penalty,
                   alpha, gamma, tau);
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
    prog.increment();
  } //end loop over lambda values
  
  double d = solver->get_d();
  
  delete solver;
  solver = NULL;
  
  return Rcpp::List::create(Named("beta")       = Rcpp::wrap(beta),
                            Named("lambda")     = lambda_tmp,
                            Named("niter")      = niter,
                            Named("d")          = d);
  
}


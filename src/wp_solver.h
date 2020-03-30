#ifndef WP_SOLVER_H
#define WP_SOLVER_H

#include <RcppEigen.h>
#include "limbs_types.h"
#include "utils.h"
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseSymMatProd.h"

class WpSolver
{
protected:
  const int nvars;                  // dimension of beta
  const int betadim;                // vector space of beta
  const int nobs;                         // number of rows
  const int nsamps;
  const int ngroups;                // number of groups for group lasso
  double power;                     // wasserstein power >=1
  
  bool intercept;                   //
  bool standardize;                 //
  
  double meanY;
  double scaleY;
  
  vector u;                       // u vector
  
  vector beta;                 // parameters to be optimized
  vector beta_prev;            // parameters from previous iteration
  vector beta_prev_irls;       // parameters from previous irls iteration
  
  double tol;                       // tolerance for convergence
  // typedef float Scalar;
  // typedef double Double;
  // typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  // typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  // typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
  // typedef Eigen::Map<const Eigen::vector> MapVec;
  // typedef Eigen::Map<const Eigen::MatrixXd> MapMatd;
  // typedef Eigen::Map<const Eigen::vector> MapVecd;
  // typedef Eigen::Map<Eigen::vectorI> MapVeci;
  // typedef const Eigen::Ref<const Eigen::MatrixXd> ConstGenericMatrix;
  // typedef const Eigen::Ref<const Eigen::vector> ConstGenericVector;
  // typedef Eigen::Ref<Eigen::vector> GenericVector;
  // typedef Eigen::Ref<Eigen::MatrixXd> GenericMatrix;
  // typedef Eigen::SparseMatrix<double> SpMat;
  // typedef Eigen::SparseVector<double> SparseVector;
  
  matrix X;                 // X matrix with scaling
  matMapConst X_orig;           // X unchanged
  SpMat X_sp;
  vecMapConst Y;                // Y vector/matrix
  vector XY;
  vector obs_weights;
  vectorI groups;            // vector of group membersihp indexes
  vectorI unique_groups;     // vector of all unique groups
  vector penalty_factor;    // penalty multiplication factors
  vector group_weights;     // group lasso penalty multiplication factors
  vector scale_factor;      // scaling factor for columns of X
  vector scale_factor_inv;  // inverse of scaling factor for columns of X
  int penalty_factor_size;    // size of penalty_factor vector
  bool selection;
  
  // Eigen::MatrixXd A;                 // A = d * I - X'X
  Eigen::SparseMatrix<double> A;
  double d;                   // d value (largest eigenvalue of X'X) for given set of weights
  bool default_group_weights; // do we need to compute default group weights?
  
  
  std::vector<std::vector<int> > grp_idx; // vector of vectors of the indexes for all members of each group
  std::string penalty;        // penalty specified
  
  double lambda;              // L1 penalty
  double lambda0;             // minimum lambda to make coefficients all zero
  double alpha;               // alpha = mixing parameter for elastic net
  double gamma;               // extra tuning parameter for mcp/scad
  double tau;                 // mixing parameter for group sparse penalties
  
  double threshval;
  int scale_len;
  
  bool found_grp_idx;
  bool is_projection;
  
  void soft_threshold(vector &res, const vector &vec, const double &penalty,
                           vector &pen_fact, double &d);
  
  
  void soft_threshold_mcp(vector &res, const vector &vec, const double &penalty,
                               vector &pen_fact, double &d, double &gamma);
  
  void soft_threshold_scad(vector &res, const vector &vec, const double &penalty,
                                vector &pen_fact, double &d, double &gamma);
  
  double soft_threshold_scad_norm(double &b, const double &pen, double &d, double &gamma);
  
  double soft_threshold_mcp_norm(double &b, const double &pen, double &d, double &gamma);
  
  void block_soft_threshold_scad(vector &res, const vector &vec, const double &penalty,
                                      vector &pen_fact, double &d,
                                      std::vector<std::vector<int> > &grp_idx,
                                      const int &ngroups, vectorI &unique_grps, vectorI &grps,
                                      double & gamma);
  
  void block_soft_threshold_mcp(vector &res, const vector &vec, const double &penalty,
                                     vector &pen_fact, double &d,
                                     std::vector<std::vector<int> > &grp_idx,
                                     const int &ngroups, vectorI &unique_grps, vectorI &grps,
                                     double & gamma);
  
  void block_soft_threshold(vector &res, const vector &vec, const double &penalty,
                                 vector &pen_fact, double &d,
                                 std::vector<std::vector<int> > &grp_idx,
                                 const int &ngroups, vectorI &unique_grps, vectorI &grps);
  
  double soft_threshold_norm (double &norm, const double &pen);
  
  void get_group_indexes();
  
  void compute_XtX_d_update_A();
  
  void compute_weights();
  
  void next_u(vector &res);
  
  void next_beta(vector &res);
  
  bool converged()
  {
    return (stopRule(beta, beta_prev, tol)); //diff
  }
  
  bool converged_irls()
  {
    return (stopRule(beta, beta_prev_irls, tol)); //diff
  }
  
public:
  WpSolver(  int s_, //diff
             const refMat  &X_,
             const refMat &Y_,
             const double &power_,
             const vectorI &groups_,
             const vectorI &unique_groups_,
             vector &group_weights_,
             vector &penalty_factor_,
             const vector &scale_factor_,
             const double tol_ = 1e-6) :
  nvars(X_.rows()*s_),
  betadim(X_.rows()), 
  nobs(Y_.rows()),
  nsamps(s_),
  ngroups(unique_groups_.size()),
  power(power_),
  obs_weights(X_.cols() * s_),
  tol(tol_),
  u(X_.rows() * s_),               // allocate space but do not set values //diff
  beta(X_.rows() * s_),            // allocate space but do not set values //diff
  beta_prev(X_.rows() * s_),       // allocate space but do not set values //diff
  beta_prev_irls(X_.rows() * s_),
                        X_orig(X_.data(), X_.rows(), X_.cols()),
                        X_sp(X_.cols(), X_.rows()),
                        X(X_.cols(), X_.rows()),
                        Y(Y_.data(), Y_.rows() * Y_.cols()),
                        XY(X_.rows() * s_),
                        groups(groups_),
                        unique_groups(unique_groups_),
                        penalty_factor(penalty_factor_),
                        group_weights(group_weights_),
                        scale_factor(scale_factor_),
                        scale_factor_inv(X_.rows()),
                        penalty_factor_size(penalty_factor_.size()),
                        default_group_weights(bool(group_weights_.size() < 1)), // compute default weights if none given
                                                                           grp_idx(unique_groups_.size())
                                                                           {}
  
  void init_wpsolve();
  
  
  double compute_lambda_zero(std::string penalty_);
  
  double get_d() { return d; }
  
  // init() is a cold start for the first lambda
  void init(double lambda_, std::string penalty_,
            double alpha_, double gamma_, double tau_);
  
  void beta_zeros() //diff
  {
    beta.setZero();
  }
  
  void init_warm(double lambda_)
  {
    lambda = lambda_;
    
  }
  
  vector get_beta();
  
  virtual double get_loss()
  {
    return 1e99;
  }
  
  void update_u()
  {
    //TypeBeta newbeta(nvars);
    next_u(u);
    //beta.swap(newbeta);
  }
  
  void update_beta()
  {
    //TypeBeta newbeta(nvars);
    next_beta(beta);
    //beta.swap(newbeta);
  }
  
  virtual void solve_param(int maxit)
  {
    
    for(int i = 0; i < maxit; ++i)
    {
      // Rcpp::Rcout << "iteration " << i << "\n";
      if(i % 1000)  Rcpp::checkUserInterrupt(); 
      
      beta_prev = beta;
      
      update_u();
      
      update_beta();
      
      if(converged()) break;
      
    }
    
  }
  
  
  virtual int solve(int maxit)
  {
    int i;
    
    for(i = 0; i < maxit; ++i)
    {
      // Rcpp::Rcout << "iteration " << i << "\n";
      if(i % 200)  Rcpp::checkUserInterrupt(); 
      
      beta_prev_irls = beta;
      
      solve_param(maxit);
      
      if(converged_irls()) break;
      
      compute_weights();
      XY = X_sp.transpose() * obs_weights.asDiagonal() * Y;
      compute_XtX_d_update_A();
      
    }
    
    
    return i + 1;
  }
  
};

#endif // WP_SOLVER

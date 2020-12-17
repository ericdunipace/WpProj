#include "SparsePosterior_types.h"
#include "utils.h"
#include "Spectra/SymEigsSolver.h"

class w2ellipse {
protected:
  // typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;
  // typedef Eigen::Matrix<int, Eigen::Dynamic, 1> vectorI;
  // 
  // typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
  // typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrixI;
  // 
  // typedef Eigen::Ref<Eigen::MatrixXd> refMat;
  // typedef Eigen::Ref<matrixI> refMatI;
  // 
  // typedef Eigen::Ref<Eigen::ArrayXi> refArrayI;
  // typedef Eigen::Ref<Eigen::ArrayXd> refArray;
  // 
  // typedef Eigen::Ref<vector> refVec;
  // 
  // typedef Eigen::Ref<const matrix> refMatConst;
  // typedef Eigen::Ref<const matrixI> refMatConstI;
  // typedef Eigen::Ref<const Eigen::ArrayXd> refArrayConst;
  // typedef Eigen::Ref<const Eigen::ArrayXi> refArrayConstI;
  // typedef Eigen::Ref<const vector> refVecConst;
  
  
  vector means; //mean vector
  vector means_old; //save last value
  vector u_means;
  vector sd; //scale matrix (lower tri)
  vector corr;
  vector Y_means;
  matrix Y_cov;
  // int p; //norm
  int d; //dimension of space
  int S; //number of samples
  int n; //observations of data
  int maxit; 
  int maxEig;
  int scale_len;
  
  Eigen::MatrixXd A;
  const matMapConst X;
  const matMapConst Y;
  const matMapConst theta;
  Eigen::MatrixXd XX;
  Eigen::MatrixXd XY_mean;
  
  Eigen::VectorXd scale_factor;
  Eigen::VectorXd scale_factor_inv;
  std::string penalty;
  
  void compute_XtX_d_update_A();
  
  void solve_means();
  void next_u_mean();
  void next_mean();
  void solve_scale();
  bool converged();
  
public:
  w2ellipse(//int p_,
            int d_, int n_,
            int S_,
            const refMatConst  &X_, 
            const refMatConst & Y_,
            const refMatConst & theta_,
            const VectorXd &scale_factor_) : 
              // p(p_),
              d(d_),
              n(n_),
              S(S_),
              X(X_.data(), X_.rows(), X_.cols()),
              Y(Y_.data(), Y_.rows(), Y_.cols()),
              theta(theta_.data(), theta_.rows(), theta_.cols()),
              XX(X_.cols(), X_.cols()),
              scale_factor(scale_factor_),
              scale_factor_inv(X_.cols()),
              corr(int( X_.cols() * (X_.cols()+1)/2.0- X_.cols() ) ){
    int corr_idx = 0;
    means = theta.rowwise().mean();
    matrix cov = covariance(theta, means);
    sd = cov.diagonal().array().sqrt();
    matrix cor = sd.asDiagonal() * cov * sd.asDiagonal();
    for(int i = 1; i < d; i ++ ){
      for (int j = 0; j < i; j ++) {
        corr(corr_idx) = cor(i,j);
        corr_idx++;
      }
    }
    Y_means = X.ldlt().solve(Y.colwise().mean());
    matrix covY = covariance(Y);
    matrix Y_L = covY.llt().matrixL();
    Y_cov = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.ldlt().solve(Y_L));
    Y_cov /= double(S_);
    XX = matrix(d, d).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X);
    matrix XY_mean_init = X * X.transpose() * Y_means;
    // scale_factor = scale_factor_;
    scale_len = scale_factor.size();
    
    if (scale_len)
    {
      scale_factor_inv = 1.0 / scale_factor.array();
      // X = scale_factor_inv.asDiagonal() * X;
      XY_mean = XY_mean_init.array().colwise() * scale_factor_inv.array();
    } else
    {
      XY_mean = XY_mean_init;
    }
    // compute XtX or XXt (depending on if n > p or not)
    // and compute A = dI - XtX (if n > p)
    compute_XtX_d_update_A();
  }
  
  vector get_means();
  matrix get_scale();
  void solve();
  
  
};

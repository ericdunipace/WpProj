#include "W2_ellipse.h"

void w2ellipse::compute_XtX_d_update_A()
{
  MatrixXd XXmat(XX.rows(), XX.cols());
  if (scale_len)
  {
    XXmat = scale_factor_inv.asDiagonal() * XX * scale_factor_inv.asDiagonal();
  } else
  {
    XXmat = XX;
  }
  Spectra::DenseSymMatProd<double> op(XXmat);
  int ncv = 4;
  if (XX.cols() < 4)
  {
    ncv = XX.cols();
  }
  
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, 1, ncv);
  
  eigs.init();
  eigs.compute(10000, 1e-10);
  vector eigenvals = eigs.eigenvalues();
  maxEig = eigenvals[0] * 1.005; // multiply by an increasing factor to be safe
  
  A = -XXmat;
  
  
  A.diagonal().array() += maxEig;
  
}

void w2ellipse::solve_means() {
  next_u_mean();
  means = u_means/maxEig;
}

void w2ellipse::next_u_mean() {
  u_means = XY_mean;
  u_means += A * means_old;
}

// void w2ellipse::threshold(std::string &penalty) {
//   if (penalty == "lasso")
//   {
//     soft_threshold(beta, u, lambda, penalty_factor, maxEig);
//   } else if (penalty == "ols")
//   {
//     beta = u / maxEig;
//   } else if (penalty == "elastic.net")
//   {
//     double denom = maxEig + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     soft_threshold(beta, u, lam, penalty_factor, denom);
//   } else if (penalty == "scad")
//   {
//     soft_threshold_scad(beta, u, lambda, penalty_factor, d, gamma);
//     
//   } else if (penalty == "scad.net")
//   {
//     double denom = d + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     if (alpha == 0)
//     {
//       lam   = 0;
//       denom = d + lambda;
//     }
//     
//     soft_threshold_scad(beta, u, lam, penalty_factor, denom, gamma);
//     
//   } else if (penalty == "mcp")
//   {
//     soft_threshold_mcp(beta, u, lambda, penalty_factor, d, gamma);
//   } else if (penalty == "mcp.net")
//   {
//     double denom = d + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     soft_threshold_mcp(beta, u, lam, penalty_factor, denom, gamma);
//     
//   } else if (penalty == "grp.lasso")
//   {
//     block_soft_threshold(beta, u, lambda, group_weights,
//                          d, grp_idx, ngroups,
//                          unique_groups, groups);
//   } else if (penalty == "grp.lasso.net")
//   {
//     double denom = d + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     block_soft_threshold(beta, u, lam, group_weights,
//                          denom, grp_idx, ngroups,
//                          unique_groups, groups);
//     
//   } else if (penalty == "grp.mcp")
//   {
//     block_soft_threshold_mcp(beta, u, lambda, group_weights,
//                              d, grp_idx, ngroups,
//                              unique_groups, groups, gamma);
//   } else if (penalty == "grp.scad")
//   {
//     block_soft_threshold_scad(beta, u, lambda, group_weights,
//                               d, grp_idx, ngroups,
//                               unique_groups, groups, gamma);
//   } else if (penalty == "grp.mcp.net")
//   {
//     double denom = d + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     
//     block_soft_threshold_mcp(beta, u, lam, group_weights,
//                              denom, grp_idx, ngroups,
//                              unique_groups, groups, gamma);
//   } else if (penalty == "grp.scad.net")
//   {
//     double denom = d + (1.0 - alpha) * lambda;
//     double lam = lambda * alpha;
//     
//     block_soft_threshold_scad(beta, u, lam, group_weights,
//                               denom, grp_idx, ngroups,
//                               unique_groups, groups, gamma);
//   } else if (penalty == "sparse.grp.lasso") {
//     double lam_grp = (1.0 - tau) * lambda;
//     double lam_l1  = tau * lambda;
//     
//     double fact = 1.0;
//     
//     // first apply soft thresholding
//     // but don't divide by d
//     soft_threshold(beta, u, lam_l1, penalty_factor, fact);
//     
//     MatrixXd beta_tmp = beta;
//     
//     // then apply block soft thresholding
//     block_soft_threshold(beta, beta_tmp, lam_grp,
//                          group_weights,
//                          d, grp_idx, ngroups,
//                          unique_groups, groups);
//   }
// } // end threshold

void w2ellipse::solve_scale() {
  means = XY_mean;
  means += A * means_old;
}

vector w2ellipse::get_means() {
  return means;
}

matrix w2ellipse::get_scale() {
  matrix scale = matrix::Zero(d,d);
  int corr_idx = 0;
  
  scale.diagonal() = sd;
  for(int i = 1; i < d; i ++ ){
    for (int j = 0; j < i; j ++) {
      scale(i,j) = corr(corr_idx) * sd(i) * sd(j);
      corr_idx++;
    }
  }
  
  return scale;
}

bool w2ellipse::converged() {
  double diff = (means_old - means).array().abs().sum();
  if(diff ==0 ){
    return 1;
  } else {
    return 0;
  }
}

void w2ellipse::solve() {
  for(int ii = 0; ii < maxit; ii++) {
    // save last values
    means_old = means;
    
    //update mean vector
    solve_means();
    
    //update scale matrix
    solve_scale();
    
    //check convergence
    if(converged()) {
      break;
    }
  }
}
#include "wp_solver.h"

void WpSolver::soft_threshold(vector &res, const vector &vec, const double &penalty,
                    vector &pen_fact, double &d)
{
  int v_size = vec.size();
  res.setZero();

  const double *ptr = vec.data();
  for(int i = 0; i < v_size; i++)
  {
    double total_pen = pen_fact(i) * penalty;

    if(ptr[i] > total_pen)
      res(i) = (ptr[i] - total_pen)/d;
    else if(ptr[i] < -total_pen)
      res(i) = (ptr[i] + total_pen)/d;
  }
}


void WpSolver::soft_threshold_mcp(vector &res, const vector &vec, const double &penalty,
                        vector &pen_fact, double &d, double &gamma)
{
  int v_size = vec.size();
  res.setZero();
  double gammad = gamma * d;
  double d_minus_gammainv = d - 1.0 / gamma;


  const double *ptr = vec.data();
  for(int i = 0; i < v_size; i++)
  {
    double total_pen = pen_fact(i) * penalty;

    if (std::abs(ptr[i]) > gammad * total_pen){
      res(i) = ptr[i]/d;
    } else if(ptr[i] > total_pen)
      res(i) = (ptr[i] - total_pen)/(d_minus_gammainv);
    else if(ptr[i] < -total_pen)
      res(i) = (ptr[i] + total_pen)/(d_minus_gammainv);

  }

}

void WpSolver::soft_threshold_scad(vector &res, const vector &vec, const double &penalty,
                         vector &pen_fact, double &d, double &gamma)
{
  int v_size = vec.size();
  res.setZero();
  double gammad = gamma * d;
  double gamma_minus1_d = (gamma - 1.0) * d;

  const double *ptr = vec.data();
  for(int i = 0; i < v_size; i++)
  {
    double total_pen = pen_fact(i) * penalty;

    if (std::abs(ptr[i]) > gammad * total_pen)
      res(i) = ptr[i]/d;
    else if (std::abs(ptr[i]) > (d + 1.0) * total_pen)
    {
      double gam_ptr = (gamma - 1.0) * ptr[i];
      double gam_pen = gamma * total_pen;
      if(gam_ptr > gam_pen)
        res(i) = (gam_ptr - gam_pen)/(gamma_minus1_d - 1.0);
      else if(gam_ptr < -gam_pen)
        res(i) = (gam_ptr + gam_pen)/(gamma_minus1_d - 1.0);
    }
    else if(ptr[i] > total_pen)
      res(i) = (ptr[i] - total_pen)/d;
    else if(ptr[i] < -total_pen)
      res(i) = (ptr[i] + total_pen)/d;

  }
}

double WpSolver::soft_threshold_scad_norm(double &b, const double &pen, double &d, double &gamma)
{
  double retval = 0.0;

  double gammad = gamma * d;
  double gamma_minus1_d = (gamma - 1.0) * d;

  if (std::abs(b) > gammad * pen)
    retval = 1.0;
  else if (std::abs(b) > (d + 1.0) * pen)
  {
    double gam_ptr = (gamma - 1.0);
    double gam_pen = gamma * pen / b;
    if(gam_ptr > gam_pen)
      retval = d * (gam_ptr - gam_pen)/(gamma_minus1_d - 1.0);
    else if(gam_ptr < -gam_pen)
      retval = d * (gam_ptr + gam_pen)/(gamma_minus1_d - 1.0);
  }
  else if(b > pen)
    retval = (1.0 - pen / b);
  else if(b < -pen)
    retval = (1.0 + pen / b);
  return retval;
}

double WpSolver::soft_threshold_mcp_norm(double &b, const double &pen, double &d, double &gamma)
{
  double retval = 0.0;

  double gammad = gamma * d;
  double d_minus_gammainv = d - 1.0 / gamma;

  if (std::abs(b) > gammad * pen)
    retval = 1.0;
  else if(b > pen)
    retval = d * (1.0 - pen / b)/(d_minus_gammainv);
  else if(b < -pen)
    retval = d * (1.0 + pen / b)/(d_minus_gammainv);

  return retval;
}

void WpSolver::block_soft_threshold_scad(vector &res, const vector &vec, const double &penalty,
                               vector &pen_fact, double &d,
                               std::vector<std::vector<int> > &grp_idx,
                               const int &ngroups, vectorI &unique_grps, vectorI &grps,
                               double & gamma)
{
  //int v_size = vec.size();
  res.setZero();

  for (int g = 0; g < ngroups; ++g)
  {
    double thresh_factor;
    std::vector<int> gr_idx = grp_idx[g];

    if (unique_grps(g) == 0) // the 0 group represents unpenalized variables
    {
      thresh_factor = 1.0;
    } else {
      double ds_norm = 0.0;
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        double val = vec(c_idx);
        ds_norm += val * val;
      }
      ds_norm = std::sqrt(ds_norm);
      // double grp_wts = sqrt(gr_idx.size());
      double grp_wts = pen_fact(g);
      //thresh_factor = std::max(0.0, 1.0 - penalty * grp_wts / (ds_norm) );
      thresh_factor = soft_threshold_scad_norm(ds_norm, penalty * grp_wts, d, gamma);
    }
    if (thresh_factor != 0.0)
    {
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        res(c_idx) = vec(c_idx) * thresh_factor / d;
      }
    }
  }
}

void WpSolver::block_soft_threshold_mcp(vector &res, const vector &vec, const double &penalty,
                              vector &pen_fact, double &d,
                              std::vector<std::vector<int> > &grp_idx,
                              const int &ngroups, vectorI &unique_grps, vectorI &grps,
                              double & gamma)
{
  //int v_size = vec.size();
  res.setZero();

  for (int g = 0; g < ngroups; ++g)
  {
    double thresh_factor;
    std::vector<int> gr_idx = grp_idx[g];

    if (unique_grps(g) == 0) // the 0 group represents unpenalized variables
    {
      thresh_factor = 1.0;
    } else {
      double ds_norm = 0.0;
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        double val = vec(c_idx);
        ds_norm += val * val;
      }
      ds_norm = std::sqrt(ds_norm);
      // double grp_wts = sqrt(gr_idx.size());
      double grp_wts = pen_fact(g);
      //thresh_factor = std::max(0.0, 1.0 - penalty * grp_wts / (ds_norm) );
      thresh_factor = soft_threshold_mcp_norm(ds_norm, penalty * grp_wts, d, gamma);
    } // fi unique grp not equal 0
    if (thresh_factor != 0.0)
    {
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        res(c_idx) = vec(c_idx) * thresh_factor / d;
      } // for looping through values in particular group
    } //fi thresh_factor != 0.0
  }// for looping through groups
}

void WpSolver::block_soft_threshold(vector &res, const vector &vec, const double &penalty,
                          vector &pen_fact, double &d,
                          std::vector<std::vector<int> > &grp_idx,
                          const int &ngroups, vectorI &unique_grps, vectorI &grps)
{
  //int v_size = vec.size();
  res.setZero();

  for (int g = 0; g < ngroups; ++g)
  {
    double thresh_factor;
    std::vector<int> gr_idx = grp_idx[g];
    /*
     for (int v = 0; v < v_size; ++v)
     {
     if (grps(v) == unique_grps(g))
     {
     gr_idx.push_back(v);
     }
     }
     */
    if (unique_grps(g) == 0)
    {
      thresh_factor = 1.0;
    } else
    {
      double ds_norm = 0.0;
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        double val = vec(c_idx);
        ds_norm += val * val;
      }
      ds_norm = std::sqrt(ds_norm);
      // double grp_wts = sqrt(gr_idx.size());
      double grp_wts = pen_fact(g);
      thresh_factor = std::max(0.0, 1.0 - penalty * grp_wts / (ds_norm) );
    }
    if (thresh_factor != 0.0)
    {
      for (std::vector<int>::size_type v = 0; v < gr_idx.size(); ++v)
      {
        int c_idx = gr_idx[v];
        res(c_idx) = vec(c_idx) * thresh_factor / d;
      }
    }
  }
}
double WpSolver::soft_threshold_norm (double &norm, const double &pen)
{
  return std::max(0.0, 1.0 - (pen / norm) );
}

void WpSolver::get_group_indexes()
{
  // if the group is any group penalty
  std::string grptxt("grp");
  if (penalty.find(grptxt) != std::string::npos)
  {
    found_grp_idx = true;
    grp_idx.reserve(ngroups);
    for (int g = 0; g < ngroups; ++g)
    {
      // find all variables in group number g
      std::vector<int> idx_tmp;
      for (int v = 0; v < betadim; ++v)
      {
        if (groups(v) == unique_groups(g))
        {
          idx_tmp.push_back(v);
        }
      }
      grp_idx[g] = idx_tmp;
    }
    // if group weights were not specified,
    // then set the group weight for each
    // group to be the sqrt of the size of the
    // group
    if (default_group_weights)
    {
      group_weights.resize(ngroups);
      for (int g = 0; g < ngroups; ++g)
      {
        if(unique_groups(g) == 0){
          group_weights(g) = 0;
        } else {
          group_weights(g) = std::sqrt(double(grp_idx[g].size()));
        }
      }
    }
  }
}

void WpSolver::init_wpsolve()
{
  scale_len = scale_factor.size();

  found_grp_idx = false;

  if (scale_len)
  {
    scale_factor_inv = 1.0 / scale_factor.array();
  }
  X = X_orig.transpose() * scale_factor.asDiagonal();

  X_sp.reserve(nobs * nvars);
  for (int j= 0; j < betadim; j ++) {
    for(int i = 0; i < nobs; i ++) {
      for(int s = 0; s < nsamps; s ++) X_sp.insert(i + nobs*s, j + betadim * s ) = X(i,j);
    }
  }
  X_sp.makeCompressed();

  // compute XtX or XXt (depending on if n > p or not)
  // and compute A = dI - XtX (if n > p)
  obs_weights = vector::Ones(nobs);
  compute_XtX_d_update_A();

  XY = X_sp.transpose() * Y;

}

void WpSolver::compute_XtX_d_update_A()
{

  SpMat XtWX = X_sp.transpose() * obs_weights.asDiagonal() * X_sp;

  // Spectra::DenseSymMatProd<double> op(XtWX);
  Spectra::SparseSymMatProd<double> op(XtWX);
  int ncv = 4;
  if (XtWX.cols() < 4)
  {
    ncv = XtWX.cols();
  }
  // Rcpp::Rcout << ncv << "\n";
  // Rcpp::stop("safety");
  // Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> > eigs(&op, 1, ncv); //object, number eig val (nev), convergence speed >= 2*nev
  Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::SparseSymMatProd<double> > eigs(&op, 1, ncv); //object, number eig val (nev), convergence speed >= 2*nev
  // Rcpp::Rcout << "what's happening??\n";
  eigs.init();
  eigs.compute(10000, 1e-10); // values are iterations and tolerance
  vector eigenvals = eigs.eigenvalues();
  d = eigenvals[0] * 1.005; // multiply by an increasing factor to be safe

  A = -XtWX;


  A.diagonal().array() += d;

}

void WpSolver::compute_weights() {
  vector res = Y;
  res -= X_sp * beta;
  obs_weights = res.array().abs().pow(power - 2.0);
  if(power == 1.0) {
    for(int n = 0; n < nvars; n ++) if( obs_weights(n) > 10000.0) obs_weights(n) = 10000.0;
  }
}

double WpSolver::compute_lambda_zero(std::string penalty_)
{ //maybe change this
  // lambda0 = XY.cwiseAbs().maxCoeff();
  int temp_size = XY.rows();
  vector temp(temp_size);

  if (!found_grp_idx) {
    penalty = penalty_;
    get_group_indexes();
  }

  if (found_grp_idx) {
    temp.resize(ngroups);
    temp.fill(0.0);


    for ( int g = 0; g < ngroups; g++ ) {
      std::vector<int> gr_idx = grp_idx[g];
      for ( int i = 0; i < gr_idx.size(); i++ ) {
        double val = XY(gr_idx[i]);
        temp(g) += val * val;
      }
    }
    temp = temp.cwiseSqrt().eval();
    temp = temp.cwiseQuotient(group_weights);
  } else {
    temp = XY.cwiseQuotient(penalty_factor).cwiseAbs().eval();
  }
  if ( penalty_factor.cwiseEqual(0.0).any() ) {
    for ( int i = 0; i < temp_size; i ++ ) {
      if ( penalty_factor(i) < Eigen::NumTraits<double>::dummy_precision() ) temp(i) = 0.0;
    }
  }
  lambda0 = temp.maxCoeff();

  return lambda0;
}


void WpSolver::init(double lambda_, std::string penalty_,
                      double alpha_, double gamma_, double tau_)
{
  beta.setZero();

  lambda = lambda_;
  penalty = penalty_;

  alpha = alpha_;
  gamma = gamma_;
  tau   = tau_;

  // get indexes of members of each group.
  // best to do just once in the beginning
  if (!found_grp_idx)
  {
    get_group_indexes();
  }
  // Rcpp::Rcout << penalty <<"\n";
  // Rcpp::Rcout << penalty <<"\n";
}


vector WpSolver::get_beta()
{
  vector res = beta;
  if (scale_len) {
    int iters = nvars/betadim;
    for(int n = 0; n < iters; n ++)
      res.segment(n,betadim).array() /= scale_factor_inv.array();
  }
  
  return res;
}

void WpSolver::next_u(vector &res)
{
  res = XY;
  res += A * beta_prev;
  // res.noalias() += A * beta_prev;
}

void WpSolver::next_beta(vector &res)
{
  if (penalty == "lasso")
  {
    soft_threshold(beta, u, lambda, penalty_factor, d);
  } else if (penalty == "ols")
  {
    beta = u / d;
  } else if (penalty == "elastic.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    soft_threshold(beta, u, lam, penalty_factor, denom);
  } else if (penalty == "scad")
  {
    soft_threshold_scad(beta, u, lambda, penalty_factor, d, gamma);
    
  } else if (penalty == "scad.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    if (alpha == 0)
    {
      lam   = 0;
      denom = d + lambda;
    }
    
    soft_threshold_scad(beta, u, lam, penalty_factor, denom, gamma);
    
  } else if (penalty == "mcp")
  {
    soft_threshold_mcp(beta, u, lambda, penalty_factor, d, gamma);
  } else if (penalty == "mcp.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    soft_threshold_mcp(beta, u, lam, penalty_factor, denom, gamma);
    
  } else if (penalty == "grp.lasso")
  {
    block_soft_threshold(beta, u, lambda, group_weights,
                         d, grp_idx, ngroups,
                         unique_groups, groups);
  } else if (penalty == "grp.lasso.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    block_soft_threshold(beta, u, lam, group_weights,
                         denom, grp_idx, ngroups,
                         unique_groups, groups);
    
  } else if (penalty == "grp.mcp")
  {
    block_soft_threshold_mcp(beta, u, lambda, group_weights,
                             d, grp_idx, ngroups,
                             unique_groups, groups, gamma);
  } else if (penalty == "grp.scad")
  {
    block_soft_threshold_scad(beta, u, lambda, group_weights,
                              d, grp_idx, ngroups,
                              unique_groups, groups, gamma);
  } else if (penalty == "grp.mcp.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    
    block_soft_threshold_mcp(beta, u, lam, group_weights,
                             denom, grp_idx, ngroups,
                             unique_groups, groups, gamma);
  } else if (penalty == "grp.scad.net")
  {
    double denom = d + (1.0 - alpha) * lambda;
    double lam = lambda * alpha;
    
    block_soft_threshold_scad(beta, u, lam, group_weights,
                              denom, grp_idx, ngroups,
                              unique_groups, groups, gamma);
  } else if (penalty == "sparse.grp.lasso") {
    double lam_grp = (1.0 - tau) * lambda;
    double lam_l1  = tau * lambda;
    
    double fact = 1.0;
    
    // first apply soft thresholding
    // but don't divide by d
    soft_threshold(beta, u, lam_l1, penalty_factor, fact);
    
    MatrixXd beta_tmp = beta;
    
    // then apply block soft thresholding
    block_soft_threshold(beta, beta_tmp, lam_grp,
                         group_weights,
                         d, grp_idx, ngroups,
                         unique_groups, groups);
  } else {
    Rcpp::stop("Penalty factor not found!");
  }
}

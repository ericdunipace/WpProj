

// [[Rcpp::export]]
Rcpp::List xtxyUpdate(const NumericMatrix & X_, const NumericMatrix & Y_,
                      const NumericMatrix & theta_,
                      const NumericVector & result_,
                      const double pseudo_obs = 0.0,
                      const int sampleIdx = 0,
                      const int D = 1,
                      const CharacterVector & method = "selection.variable") {
  const int S = theta_.rows();
  const int P = X_.cols();
  const int N = X_.rows();
  const double mix_wt = pseudo_obs/(double(N) + pseudo_obs);
  
  const matMap Xt(as<matMap >( X_ ));
  const matMap Yt(as<matMap >( Y_ ));
  const matMap thetat(as<matMap>( theta_ ));
  const matMap result(as<matMap>( result_ ));
  
  const matrix X = Xt.transpose();
  matrix Y = Yt.transpose();
  matrix theta = thetat.transpose();
  matrix mu(S,N);
  matrix theta_norm(P,1);
  
  matrixI idx_mu(S,N);
  matrix xty = matrix::Zero(P,1);
  matrix xtx = matrix::Zero(P,P);
  
  if ( method(0) != "projection" ) sort_matrix(Y);
  if ( method(0) == "projection") {
    // theta.resize(P,S);
    theta_norm.resize(P,S);
    xty.resize(P,S);
    xty = X * Yt * (1.0 - mix_wt)/double(N);
  } // fi projection
  if ( method(0) == "location.scale") {
    Rcpp::stop("Location/scale method not supported");
    if (result.rows() != 2*P) {
      // if(result.rows() == P) {
      //   Rcpp::warning("Doubling input result vector to match method");
      //   vector temp_result = result;
      //   result.resize(2*P,1);
      //   result.block(0,0,P,1) = temp_result;
      //   result.block(P,0,P,1) = temp_result;
      // } else {
      Rcpp::stop("Result vector must have dimension twice that of theta");
      // }
    }
    matrix means(P,S);
    
    {
      vector meanvec = theta.rowwise().mean();
      for (int s = 0; s < S ; s++) means.col(s) = meanvec;
    }
    
    matrix c_theta = theta - means;
    
    theta.resize(2*P,S);
    theta_norm.resize(2*P,1);
    xty.resize(2*P,1);
    
    theta.block(0,0,P,S) = c_theta;
    theta.block(P,0,P,S) = means;
    
  } // fi loc.scale
  xtxy_update(X, Y, theta,
              theta_norm, result,
              mu, S,  N, D,
              mix_wt,
              sampleIdx,
              xtx, xty, idx_mu,
              method);
  
  return Rcpp::List::create(Rcpp::Named("XtX") = Rcpp::wrap(xtx),
                            Rcpp::Named("XtY") = Rcpp::wrap(xty));
}

// void xty_update(const refMatConst & X, 
//                 const refMatConst & sorted_Y, //y needs to be pre-sorted for scale/loc.scale
//                 const refMatConst & theta, 
//                 matrix & theta_norm,
//                 const refMatConst & result,
//                 matrix & mu,
//                 const int S, 
//                 const int N,
//                 const double mix_wt,
//                 matrix & xty, 
//                 matrixI & idx_mu, 
//                 const Rcpp::CharacterVector & method) {
//   
//   if( (method(0) == "scale") || (method(0) == "selection.variable") ) {
//     mu_update(X, result, theta, mu, method);
//     xty_update_scale(X, sorted_Y, theta, theta_norm, result, mu, S, N, mix_wt, xty, idx_mu);
//   // } else if (method(0) == "location.scale") {
//   //   mu_update(X, result, theta, mu, method);
//   //   xty_update_loc_scale(X, sorted_Y, theta,
//   //                        theta_norm, result, mu, S, N, D, wt, sampleIdx, xty, idx_mu);
//   // } else if (method(0) == "projection") {
//     // mu_update(X, result, theta, mu, method); //only need if sorting!!
//     // xty_update_projection(X, sorted_Y, theta, theta_norm, result, //result=theta_perp
//     //                        mu, S, N, mix_wt, xty, idx_mu);
//   } else {
//     Rcpp::stop("Method not found in update xty!");
//   }
//   
// }




// void xty_update_eff(const refMatConst & X, const refMatConst & sorted_Y, //y needs to be pre-sorted
//                 const refMatConst & theta,
//                 const refMatConst &mu,
//                 const int S, const int P, const int N,
//                 const double Div_N_d,
//                 vector & xty, const matrixI & idx_mu, const matrixI & old_idx_mu,
//                 const vectorI & idx_col,
//                 matrix & xty_temp) {
// 
//   // matrix temp(P,S);
//   vector temp_Y(S);
// 
// 
//   for ( int n = 0; n < N; n++ ) {
//     if(n % 10 == 0) Rcpp::checkUserInterrupt();
//     
//     int curcol = idx_col(n);
//     if(curcol == -1) break;
//     
//     matrix::ColXpr tempcol = xty_temp.col( curcol );
//     xty.noalias() -= tempcol;
// 
//     //fill vector based on order of predicted mean and sorted Y
//     rel_sorted_1( idx_mu.col(curcol), temp_Y, sorted_Y.col(curcol) );
//     
//     //calculate XtY
//     tempcol = (theta.array().colwise() * X.col(n).array()).matrix() * Div_N_d * temp_Y;
//     xty.noalias() +=  tempcol;
//   }
// }
// void xty_update_loc_scale(const refMatConst & X, const refMatConst & sorted_Y, //y needs to be pre-sorted
//                           const refMatConst & theta,
//                           matrix & theta_norm,
//                           const refMatConst & result,
//                           const refMatConst &mu,
//                           const int S, const int N,
//                           const double wt,
//                           matrix & xty, matrixI & idx_mu) {
//   
//   vector temp_Y(S);
//   int P = theta.rows()/2;
//   matrix c_theta = theta.block(0,0,P,S);
//   matrix theta_mean = theta.block(P,0,P,S);
//   double dataWt = (1.0 - wt)/double( N * S );
//   double postWt = wt/double(S);
//   
//   sort_indexes_Eigenmat(mu, idx_mu);
// 
//   if(wt != 0) {
//     //recreate theta and get sparse theta for transport
//     matrix sparse_theta_full = result.asDiagonal() * theta;
//     matrix sparse_theta = sparse_theta_full.block(0,0,P,S) + sparse_theta_full.block(P,0,P,S);
//     matrix temp_theta(2 * P, 1);
//     // matrix recombine_theta = c_theta + theta_mean;
//     // 
//     // //data for transport
//     // matrixI assign_mat_theta = matrixI::Zero(S,S);
//     // matrixI basis_mat_theta = matrixI::Zero(S,S);
//     // matrix dist_theta = matrix::Zero(S,S);
//     // matrixI idx_theta(S,S);
//     // vector mass_theta(S);
//     // 
//     // //optimal transport for posterior
//     // transport(recombine_theta, sparse_theta, 2.0, 2.0,
//             // idx_theta, mass_theta, "shortsimplex")
//     // 
//     // for(int s = 0; s < S; s++) {
//     //   temp_theta.col(idx_theta(s,1)) = theta.col(idx_theta(s,0));
//     // }
//     temp_theta = theta.rowwise().squaredNorm() ;
//     theta_norm = temp_theta * postWt;
//     xty = theta_norm;
//   } else {
//     xty.fill(0.0);
//   }
//   
//   for ( int n = 0; n < N; n++ ) {
//     if(n % 10 == 0) Rcpp::checkUserInterrupt();
// 
//     //fill vector based on order of predicted mean and sorted Y
//     rel_sorted_1(idx_mu.col(n), temp_Y, sorted_Y.col(n) );
// 
//     //calculate XtY
//     vector xcol = X.col(n);
//     xty.block(0,0,P,1).noalias() +=  (c_theta.array().colwise() * xcol.array()).matrix() * temp_Y * dataWt;
//     xty.block(P,0,P,1).noalias() +=  (theta_mean.array().colwise() * xcol.array()).matrix() * temp_Y * dataWt;
//   }
//   
// }

// void xty_update_projection(const refMatConst & X, const refMatConst & Y, //y needs to be pre-sorted
//                       const refMatConst & theta,
//                       matrix & theta_norm,
//                       const refMatConst & theta_perp,
//                       const refMatConst &mu,
//                       const int S, const int N,
//                       const double mix_wt,
//                       matrix & xty, const matrixI & idx_mu) {
//   //weights
//   double dataWt = (1.0 - mix_wt)/double(N);
//   double postWt = mix_wt/double(S);
//   
//   //sorting variables for mean
//   // matrixI assign_mat_mu = matrixI::Zero(S,S);
//   // matrixI basis_mat_mu = matrixI::Zero(S,S);
//   // matrix dist_mu = matrix::Zero(S,S);
//   // vectorI idx = vectorI::Zero(S);
//   
//   //xty calculation for mean
//   // matrix xty_temp = X * Y.adjoint() * dataWt;
//   // xty.setZero();
//   
//   //update component for penalty for straying from posterior
//   if (mix_wt != 0.0) {
//     matrix xty_temp = X * Y.adjoint() * dataWt;
//     vectorI idx_theta = vectorI::Zero(S);
//     
//     //short simplex only good if have a distribution with more than one dimension
//     if (theta.rows() > 1){
//       // matrixI assign_mat_theta = matrixI::Zero(S,S);
//       // matrixI basis_mat_theta = matrixI::Zero(S,S);
//       // matrix dist_theta = matrix::Zero(S,S);
//       
//       //optimal transport for posterior
//       matrixI idx_theta(S*S,2);
//       vector mass(S);
//       
//       transport(theta, theta_perp, 2.0, 2.0,
//                 idx_theta, mass, "shortsimplex");
//     } else { //sorts the univariate distributions
//       const vector theta_sort = theta;
//       const vector theta_perp_sort = theta_perp;
//       const std::vector<size_t> t_idx = sort_indexes(theta);
//       const std::vector<size_t> tp_idx = sort_indexes(theta_perp);
//       for ( int s=0; s < S; s++ ) idx_theta(tp_idx[s]) = t_idx[s];
//     }
//     for (int s = 0; s < S; s++) theta_norm.col(s) = theta.col(idx_theta(s)) * postWt;
//     xty = theta_norm;
//     xty += xty_temp;
//   }
//   // vectorI idx_ord = vectorI::LinSpaced(S,0,S-1);
//   
//   // optimal transport for mu. may not be necessary for projection method
//   //remove mu_update method if not sorting!
//   // if (N > 1) {
//   //   idx = transport(Y.transpose(), mu.transpose()); //not latest function
//   // } else {
//   //   const vector mu_sort = mu;
//   //   const vector Y_sort = Y;
//   //   const std::vector<size_t> mu_idx = sort_indexes(mu_sort);
//   //   const std::vector<size_t> Y_idx = sort_indexes(Y_sort);
//   //   for ( int s=0; s < S; s++ ) idx(mu_idx[s]) = Y_idx[s];
//   // }
//   // Rcout << idx << std::endl;
//   // Rcout << idx_ord.cwiseEqual(idx).count() <<", if true " << S << std::endl;
//   
//   // just need to permute the columns of xty based on infimum found
//   // double W2 = 0.0;
//   // for ( int s=0; s < S; s++ ) {
//   //   xty.col(s) += xty_temp.col(idx(s));
//   //   W2 += ( Y.row(idx(s)) - mu.row(s) ).squaredNorm() / double(S);
//   // }
//   // Rcpp::Rcout << "\nW2^2: " << W2 <<"\n";
//   // Rcpp::Rcout << mu(0,0) << "\n";
// }

void xtxy_update_projection(const refMatConst & X, const refMatConst & Y, //y needs to be pre-sorted
                            const refMatConst & theta,
                            matrix & theta_norm,
                            const refMatConst & theta_perp,
                            const refMatConst &mu,
                            const int S, const int N, const int D,
                            const double mix_wt,
                            const int sampleIdx,
                            matrix & xtx, 
                            matrix & xty,
                            const matrixI & idx_mu) {
  int P = xtx.cols();
  //weights
  double dataWt = (1.0 - mix_wt)/double(N);
  // double postWt = pseudo_obs/double(S);
  vector resid = (Y.row(sampleIdx) - mu.row(sampleIdx)).transpose();
  vector wt = resid.array().pow(double(D) - 2.0);
  for (int n = 0 ; n< N; n++) {
    if(wt(n) == 0.0 ) wt(n) = Eigen::NumTraits<double>::epsilon();
  }
  if (D == 1) {
    for(int i = 0; i < wt.size(); i++) {
      if(wt(i) == 0.0) wt(i) = Eigen::NumTraits<double>::epsilon();
    }
  }
  vector wt_half = wt.cwiseSqrt();
  xtx = matrix( P, P ).setZero().selfadjointView<Lower>().rankUpdate(X * wt_half.asDiagonal(), dataWt);
  xty.col(sampleIdx) = X * wt.asDiagonal() * Y.row(sampleIdx).transpose();
  
}

// void suff_stat_loc_scale(const refMatConst & X, refMat Y,
//                          refMat theta,
//                          const bool not_same,
//                          const int S, const int P, const int N, const int D,
//                          const double pseudo_obs,
//                          matrix & xtx_dens, matrix & xty, matrix & theta_norm) {
//   double dataWt = (1.0-pseudo_obs)/double(N*S);
//   double postWt = pseudo_obs/double(S);
//   matrix xtx = matrix::Zero(2 * P, 2 * P);
//   matrix temp(2 * P, S);
//   vector temp_Y(S);
//   vector mu(S);
//   std::vector<size_t> idx_mu(S);
//   
//   if (xtx_dens.size() != (2 * P * 2 * P) ){
//     xtx_dens.resize(2 * P, 2 * P);
//   }
//   if (xty.size() != (2 * P)) {
//     xty.resize(2 * P, 1);
//   }
//   theta_norm.resize(2 * P, 1);
//   
//   matrix means(P,S);
//   
//   {
//     vector meanvec = theta.rowwise().mean();
//     for (int s = 0; s < S ; s++) means.col(s) = meanvec;
//   }
//   
//   matrix c_theta = theta - means;
//   
//   {
//     matrix c_theta_prod = c_theta.rowwise().squaredNorm() * postWt;
//     matrix m_theta_prod = means.col(0);
//     theta_norm.block(0,0,P,1) = c_theta_prod;//c_theta_prod.cwiseSqrt(); // if parameterize directly as sd
//     theta_norm.block(P,0,P,1) = m_theta_prod.cwiseAbs2() * wt;
//   }
//   
//   theta_norm.array() *= postWt;
//   xty = theta_norm;
//   
//   for ( int n = 0; n < N; n++ ) {
//     if( n % 10 == 0 ) Rcpp::checkUserInterrupt();
//     vector xcol = X.col(n);
//     
//     //Block of size (p,q), starting at (i,j)	matrix.block(i,j,p,q); 
//     //https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
//     temp.block(0, 0, P, S) = c_theta.array().colwise() * xcol.array();
//     temp.block(P, 0, P, S) =   means.array().colwise() * xcol.array();
//     temp_Y = Y.col(n);
//     
//     if ( not_same ) {
//       mu = temp.colwise().sum();
//       sort_indexes(mu, idx_mu);
//       rel_sort(idx_mu, temp_Y);
//     }
//     xtx.selfadjointView<Eigen::Lower>().rankUpdate(temp , dataWt);
//     // if(n == 0) {
//     //   writeToCSVfile("c_theta.txt",c_theta);
//     //   writeToCSVfile("m_theta.txt", means);
//     //   writeToCSVfile("temp.txt", temp);
//     //   writeToCSVfile("theta.txt", theta);
//     // }
//     xty.noalias() +=  temp * temp_Y * dataWt;
//   }
//   xtx_dens = xtx.selfadjointView<Eigen::Lower>();
//   // writeToCSVfile("xtx1.txt", xtx_dens);
//   // matrix diags = theta_norm.asDiagonal();
//   // Rcout << diags;
//   // writeToCSVfile("theta_norm.txt", theta_norm);
//   xtx_dens.diagonal() += theta_norm;
//   // writeToCSVfile("xtx2.txt", xtx_dens);
//   theta.resize(2*P, S);
//   theta.block(0,0,P,S) = c_theta;
//   theta.block(P,0,P,S) = means;
//   sort_matrix(Y);
// }

// void suff_stat_projection(const refMatConst & X, refMat Y,
//                           refMat theta,
//                           const bool not_same,
//                           const int S, const int P, const int N,
//                           const double mix_wt,
//                           matrix & xtx, matrix & xty, matrix & theta_norm) {
//   if(xty.cols() != S) xty.resize(P,S);
//   
//   double dataWt = ( 1.0 - mix_wt ) / double( N );
//   double postWt = mix_wt/double(S);
//   
//   xty = X * Y.adjoint() * dataWt;
//   xtx = matrix( P, P ).setZero().selfadjointView<Lower>().rankUpdate(X , dataWt);
//   if( mix_wt != 0.0) {
//     theta_norm = theta * postWt;
//     xty += theta_norm;
//     xtx += matrix::Identity(P,P) * postWt;
//   }
//   
// }

// void writeToCSVfile(std::string fn, matrix matrix)
// {
//   std::ofstream file(fn.c_str());
//   if (file.is_open())
//   {
//     Eigen::IOFormat full(Eigen::FullPrecision,0," ", "\n","","","","");
//     file << matrix.format(full) << '\n';
//     //file << "m" << '\n' <<  colm(matrix) << '\n';
//   }
//   file.close();
// }
// 
// void writeToCSVfile(std::string fn, vectorI matrix)
// {
//   std::ofstream file(fn.c_str());
//   if (file.is_open())
//   {
//     // Eigen::IOFormat full(Eigen::FullPrecision,0," ", "\n","","","","");
//     file << matrix << '\n';
//     //file << "m" << '\n' <<  colm(matrix) << '\n';
//   }
//   file.close();
// }

// void suff_stat_scale(const refMatConst & X, refMat Y,
//                      const refMatConst & theta,
//                      const bool not_same,
//                      const int S, const int P, const int N,
//                      const double mix_wt,
//                      matrix & xtx_dens, matrix & xty, matrix & theta_norm) {
//   double dataWt = (1.0 - mix_wt)/double(N * S);
//   double postWt = mix_wt/double(S);
//   matrix xtx = matrix::Zero(P,P);
//   matrix temp(P, S);
//   vector temp_Y(S);
//   vector mu(S);
//   std::vector<size_t> idx_mu(S);
//   
//   theta_norm = theta.rowwise().squaredNorm() * postWt;
//   
//   if(mix_wt != 0.0 ) {
//     xty = theta_norm;
//   } else {
//     xty.setZero();
//   }
//   
//   for ( int n = 0; n < N; n++ ) {
//     if(n % 10 == 0) Rcpp::checkUserInterrupt();
//     
//     temp = theta.array().colwise() * X.col(n).array();
//     temp_Y = Y.col(n);
//     
//     if ( not_same ) {
//       mu = temp.colwise().sum();
//       sort_indexes(mu, idx_mu);
//       rel_sort(idx_mu, temp_Y);
//     }
//       xtx.selfadjointView<Eigen::Lower>().rankUpdate(temp , dataWt);
//       xty.noalias() +=  temp * temp_Y * dataWt;
//   }
//   xtx_dens = xtx.selfadjointView<Eigen::Lower>();
//   
//   if(mix_wt != 0.0 ) {
//     xtx_dens += theta_norm.asDiagonal();
//   }
//   sort_matrix(Y);
// }










void xtxy_update_scale(const refMatConst & X, const refMatConst & Y, //y needs to be pre-sorted
                       const refMatConst & theta,
                       matrix & theta_norm,
                       const refMatConst & result,
                       const refMatConst &mu,
                       const int S, const int N, const int D,
                       const double mix_wt,
                       const int sampleIdx,
                       matrix & xtx, 
                       matrix & xty,
                       matrixI & idx_mu) {
  double dataWt = (1.0 - mix_wt)/double(N*S);
  double postWt = mix_wt/double(S);
  bool use_weights = (D != 2);
  sort_indexes_bycol_Eigenmat(mu, idx_mu);
  
  int P = theta.rows();
  vector temp_Y(S);
  matrix temp(P, S);
  // xtx.fill(0.0);
  matrix xtx_temp = matrix::Zero(P,P);
  
  for ( int n = 0; n < N; n++ ) {
    if(n % 10 == 0) Rcpp::checkUserInterrupt();
    
    temp = theta.array().colwise() * X.col(n).array();
    
    rel_sorted_1( idx_mu.col(n), temp_Y, Y.col(n) );
    
    vector resid = temp_Y - mu.col(n);
    vector wt = resid.cwiseAbs().array().pow( double(D) - 2.0 );
    vector wt_half = wt.cwiseSqrt();
    
    if(D == 1) {
      for(int w = 0; w < wt.size(); w++) {
        if(wt(w) == 0.0) wt(w) = Eigen::NumTraits<double>::epsilon();
      }
    }
    // Rcpp::Rcout << wt(0) << std::endl;
    if(use_weights) {
      xtx_temp.selfadjointView<Eigen::Lower>().rankUpdate( temp * wt_half.asDiagonal(), dataWt );
      xty.noalias() +=  temp * wt.asDiagonal() * temp_Y * dataWt;
    } else {
      xtx_temp.selfadjointView<Eigen::Lower>().rankUpdate( temp , dataWt );
      xty.noalias() +=  temp * temp_Y * dataWt;
    }
  }
  xtx = xtx_temp.selfadjointView<Eigen::Lower>();
  
}





void xtxy_update(const refMatConst & X, 
                 const refMatConst & sorted_Y, //y needs to be pre-sorted for scale/loc.scale
                 const refMatConst & theta, 
                 matrix & theta_norm,
                 const refMatConst & result,
                 matrix & mu,
                 const int S, 
                 const int N,
                 const int D,
                 const double mix_wt,
                 const int sampleIdx,
                 matrix & xtx, 
                 matrix & xty,
                 matrixI & idx_mu, 
                 const Rcpp::CharacterVector & method) {
  
  if( (method(0) == "scale") || (method(0) == "selection.variable") ) {
    mu_update(X, result, theta, mu, method);
    xtxy_update_scale(X, sorted_Y, theta, theta_norm, result, mu, S, N, D, mix_wt, sampleIdx, xtx, xty, idx_mu);
    // } else if (method(0) == "location.scale") {
    // mu_update(X, result, theta, mu, method);
    // xtx_update_loc_scale(X, sorted_Y, theta,
    // theta_norm, result, mu, S, N, D, pseudo_obs, sampleIdx, xtx, idx_mu);
  } else if (method(0) == "projection") {
    // mu_update(X, result, theta, mu, method); //only need if sorting!!
    xtxy_update_projection(X, sorted_Y, theta, theta_norm, result, //result=theta_perp
                           mu, S, N, D, mix_wt, sampleIdx, xtx, xty, idx_mu);
  } else {
    Rcpp::stop("Method not found in update xtx!");
  }
  
}



void xty_update_scale(const refMatConst & X, const refMatConst & sorted_Y, //y needs to be pre-sorted
                      const refMatConst & theta,
                      matrix & theta_norm,
                      const refMatConst & result,
                      const refMatConst &mu,
                      const int S, const int N,
                      const double mix_wt,
                      matrix & xty, matrixI & idx_mu) {
  double dataWt = (1.0 - mix_wt)/double(N*S);
  double postWt = mix_wt/double(S);
  sort_indexes_bycol_Eigenmat(mu, idx_mu);
  
  vector temp_Y(S);
  int P = theta.rows();
  
  if(mix_wt != 0) {
    matrix sparse_theta = result.asDiagonal() * theta;
    matrix temp_theta(P,S);
    
    // matrixI assign_mat_theta = matrixI::Zero(S,S);
    // matrixI basis_mat_theta = matrixI::Zero(S,S);
    // matrix dist_theta = matrix::Zero(S,S);
    matrixI idx_theta(S*S, 2);
    vector mass(S*S);
    bool a_sort = false;
    
    //optimal transport for posterior
    // transport(const refMatConst & A, const refMatConst & B, const double p, const double ground_p,
    //           refMat idx, refVec mass, const std::string & method)
    transport(theta, sparse_theta, 2.0, 2.0,
              idx_theta, mass, "shortsimplex", a_sort);
    
    //generate sorted theta_matrix
    for(int s = 0; s < S; s++) temp_theta.col(idx_theta(s,1)) = theta.col(idx_theta(s,0));
    
    // vectorI idx_ord = vectorI::LinSpaced(S,0,S-1);
    // Rcout << idx_theta.cwiseEqual(idx_ord).count() << std::endl;
    
    //do products of sorts
    temp_theta.array() *= theta.array() * postWt;
    
    //regenerate theta_norm
    theta_norm = temp_theta.rowwise().sum();
    xty = theta_norm;
  } else {
    xty.fill(0.0);
  }
  
  for ( int n = 0; n < N; n++ ) {
    if(n % 10 == 0) Rcpp::checkUserInterrupt();
    
    //fill vector based on order of predicted mean and sorted Y
    rel_sorted_1(idx_mu.col(n), temp_Y, sorted_Y.col(n) );
    
    //update xty
    matrix temp_prod = theta.array().colwise() * X.col(n).array();
    xty.noalias() +=  temp_prod * temp_Y * dataWt;
    
  }
  
}







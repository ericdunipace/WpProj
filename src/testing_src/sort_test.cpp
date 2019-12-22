#include <RcppEigen.h>

typedef Eigen::VectorXd vector;
typedef Eigen::VectorXi vectorI;

typedef Eigen::MatrixXd matrix;
typedef Eigen::MatrixXi matrixI;

typedef Eigen::Ref<matrix> refMat;
typedef Eigen::Ref<matrixI> refMatI;

typedef Eigen::Ref<Eigen::ArrayXi> refArrayI;
typedef Eigen::Ref<Eigen::ArrayXd> refArray;

typedef Eigen::Ref<vector> refVec;

typedef Eigen::Ref<const matrix> refMatConst;
typedef Eigen::Ref<const matrixI> refMatConstI;
typedef Eigen::Ref<const Eigen::ArrayXd> refArrayConst;
typedef Eigen::Ref<const Eigen::ArrayXi> refArrayConstI;
typedef Eigen::Ref<const vector> refVecConst;

typedef matrix::ColXpr ColXpr;
typedef matrixI::ColXpr ColXprI;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMat;
typedef Eigen::LLT<matrix> llt;
typedef Eigen::LDLT<matrix> ldlt;

typedef Eigen::Map<matrix> matMap;
typedef Eigen::Map<const matrix> matMapConst;
typedef Eigen::Map<rowMat> rowMatMap;

typedef Eigen::Map<Eigen::VectorXd> vecMap;
typedef Eigen::Map<const vector> vecMapConst;
typedef Eigen::Map<Eigen::VectorXi> vecMapI;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

bool compare (double a, double b) {
  return a < b;
}

std::vector<size_t> sort_indexes( matrix::ConstColXpr & v) {
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
  return idx;
}

void sort_indexes(const refVecConst & v, vectorI & idx) {
  int P = idx.size();
  // sort indexes based on comparing values in v
  std::sort(idx.data(), idx.data() + P,
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

std::vector<size_t> sort_indexes( const vector & v) {
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
  return idx;
}

void sort_indexes(const matrix::ConstColXpr & v, std::vector<size_t> & idx) {
  
  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes_col(const matrix::ColXpr & v, std::vector<size_t> & idx) {
  
  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes( vector & v, std::vector<size_t> & idx) {
  
  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes_Eigenmat(const refMatConst & v, matrixI & idx) { //checked
  
  int N = v.cols();
  int P = v.rows();
  
  for(int n = 0; n < N;  n++) {
    // initialize original index locations
    idx.col(n) = Eigen::VectorXi::LinSpaced(Eigen::Sequential,P,0,P-1); //fills with increasing values
    
    // sort indexes based on comparing values in v
    std::sort(idx.col(n).data(), idx.col(n).data() + P,
              [&v,n](size_t i1, size_t i2) {return v(i1,n) < v(i2,n);});
    // for(int i =0; i < idx.size(); i++) idx(i,n) = idx_temp[i];
  }
  
}

void sort_matrix(refMat v) {
  
  int N = v.cols();
  int P = v.rows();
  
  // for(auto n : v.colwise()) {
  //   std::sort(n.begin(), n.end())
  // } //available in future eigen.
  
  for(int n = 0; n < N;  n++) {
    // sort indexes based on comparing values in v
    std::sort(v.col(n).data(), v.col(n).data() + P);
  }
  
}

void rel_sort(const std::vector<size_t> & idx, vector & y) {
  vector temp_sort = y;
  std::sort(temp_sort.data(), temp_sort.data() + temp_sort.size());
  for(int i = 0; i < y.size(); i ++) y(idx[i]) = temp_sort(i);
}

// void rel_sort(const matrixI::ColXpr & idx, vector & y) {
//   vector temp_sort = y;
//   std::sort(temp_sort.data(), temp_sort.data() + temp_sort.size());
//   for(int i = 0; i < y.size(); i ++) y(idx(i)) = temp_sort(i);
// }

void rel_sort(const refArrayConstI & idx, vector & y) {
  vector temp_sort = y;
  if(temp_sort.size() != idx.size()) Rcpp::stop("Index and vector size do not match");
  std::sort(temp_sort.data(), temp_sort.data() + temp_sort.size());
  for(int i = 0; i < temp_sort.size(); i ++) y(idx(i)) = temp_sort(i);
}

void rel_sorted_1(const refArrayConstI&  idx, 
                  vector & y, const refArrayConst& yorig) {
  for(int i = 0; i < y.size(); i ++) y(idx(i)) = yorig(i);
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix sort_check_eig(Rcpp::NumericMatrix x)//, Rcpp::CharacterVector & method) 
{
  int N = x.rows();
  int P = x.cols();
  matrixI idx(N,P);
  const matMap X(as<matMap >( x ));
  sort_indexes_Eigenmat(X, idx);
  
  return Rcpp:: wrap(idx);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sort_check_rel(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y)//, Rcpp::CharacterVector & method) 
{
  int N = x.rows();
  int P = x.cols();
  const matMap X(as<matMap >( x ));
  const matMap Y(as<matMap >( y ));
  matrix out(y.rows(), y.cols());
  matrixI idx(N,P);
  sort_indexes_Eigenmat(X, idx);
  // std::vector<size_t> idx_best(N);
  
  for(int i = 0; i < P; i ++){
    vector temp = Y.col(i);
    rel_sort(idx.col(i), temp);
    out.col(i) = temp;
  }
  return Rcpp:: wrap(out);

}

// [[Rcpp::export]]
Rcpp::NumericMatrix sort_check_rel2(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y)//, Rcpp::CharacterVector & method) 
{
  int N = x.rows();
  int P = x.cols();
  const matMap X(as<matMap >( x ));
  const matMap Y(as<matMap >( y ));
  matrix out(y.rows(), y.cols());
  matrixI idx(N,P);
  sort_indexes_Eigenmat(X, idx);
  vector temp(N);
  matrix YY = Y;
  sort_matrix(YY);
  // std::vector<size_t> idx_best(N);
  
  for(int i = 0; i < P; i ++){
    rel_sorted_1(idx.col(i), temp, YY.col(i));
    out.col(i) = temp;
  }
  return Rcpp:: wrap(out);
  
}


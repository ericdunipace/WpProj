#ifndef PACKAGENAME_TYPEDEFS_H
#define PACKAGENAME_TYPEDEFS_H
#include <RcppEigen.h>

typedef Eigen::VectorXd vectorxd;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  1> vectorLD;
typedef Eigen::VectorXi vectorI;

typedef Eigen::MatrixXd matrix;
typedef Eigen::MatrixXi matrixI;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  Eigen::Dynamic> matrixLD;

typedef Eigen::Ref<matrix> refMat;
typedef Eigen::Ref<matrixI> refMatI;

typedef Eigen::Ref<Eigen::ArrayXi> refArrayI;
typedef Eigen::Ref<Eigen::ArrayXd> refArray;

typedef Eigen::Ref<vectorxd> refVec;
typedef Eigen::Ref<vectorI> refVecI;

typedef Eigen::Ref<const matrix> refMatConst;
typedef const Eigen::Ref<const matrix> constRefMatConst;
typedef Eigen::Ref<const matrixI> refMatConstI;
typedef Eigen::Ref<const Eigen::ArrayXd> refArrayConst;
typedef Eigen::Ref<const Eigen::ArrayXi> refArrayConstI;

typedef Eigen::Ref<const vectorxd> refVecConst;
typedef Eigen::Ref<const vectorI> refVecConstI;

typedef matrix::ColXpr ColXpr;
typedef matrixI::ColXpr ColXprI;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMat;
typedef Eigen::LLT<matrix> llt;
typedef Eigen::LDLT<matrix> ldlt;

typedef Eigen::Map<matrix> matMap;
typedef Eigen::Map<const matrix> matMapConst;
typedef Eigen::Map<Eigen::Matrix<long double, Eigen::Dynamic,  Eigen::Dynamic>> matMapLD;
typedef Eigen::Map<rowMat> rowMatMap;

typedef Eigen::Map<Eigen::VectorXd> vecMap;
typedef Eigen::Map<vectorLD> vecMapLD;
typedef Eigen::Map<const vectorxd> vecMapConst;
typedef Eigen::Map<Eigen::VectorXi> vecMapI;
typedef Eigen::Map<const Eigen::VectorXi> vecMapConstI;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;

#endif

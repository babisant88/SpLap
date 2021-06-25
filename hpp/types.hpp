/*
 * types.hpp
 *
 *  Created on: Nov 30, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_TYPES_HPP_
#define MYHPPFILES_TYPES_HPP_

//#define EIGEN_USE_MKL_ALL

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>

using namespace Eigen;

enum _bench_type_ {RC = 1, RLC = 2};

typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
typedef Matrix<double, Dynamic, 1> VectorXd;
typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

#endif /* MYHPPFILES_TYPES_HPP_ */

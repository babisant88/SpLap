/*
 * pinv.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_PINV_HPP_
#define MYHPPFILES_PINV_HPP_

#include "types.hpp"

template<typename _Matrix_Type_>
_Matrix_Type_ pinv(const _Matrix_Type_& a)
{
	double epsilon = std::numeric_limits<double>::epsilon();

	Eigen::BDCSVD < _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

#endif /* MYHPPFILES_PINV_HPP_ */

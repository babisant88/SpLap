/*
 * Sparsify.hpp
 *
 *  Created on: Jun 5, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_SPARSIFY_HPP_
#define MYHPPFILES_SPARSIFY_HPP_

#include "types.hpp"

SpMat Sparsify_top(MatrixXd &, int);
double* Sparsify(MatrixXd &, int, graph*&, bool);

#endif /* MYHPPFILES_SPARSIFY_HPP_ */

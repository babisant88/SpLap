/*
 * aux.hpp
 *
 *  Created on: Dec 9, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_AUX_HPP_
#define MYHPPFILES_AUX_HPP_

#include "types.hpp"

void printMatrixXd(MatrixXd &);
void printVectorXd(VectorXd &);
void printSpMat(SpMat &);
void printVectorXdToFile(VectorXd&, std::string);
void printMatrixXdToFile(MatrixXd&, std::string);
void printIntArrToFile(int* X, int n, std::string filename);

void find(SpMat&, double*, int*);

#endif /* MYHPPFILES_AUX_HPP_ */

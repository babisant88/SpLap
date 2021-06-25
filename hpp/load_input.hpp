/*
 * load_input.hpp
 *
 *  Created on: Jun 12, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_LOAD_INPUT_HPP_
#define MYHPPFILES_LOAD_INPUT_HPP_

#include <string>
#include "types.hpp"

using namespace std;

void load_dense_matrix(char*, MatrixXd&);
void load_input(char*[], MatrixXd&);

#endif /* MYHPPFILES_LOAD_INPUT_HPP_ */

/*
 * hist.cpp
 *
 *  Created on: Jun 7, 2018
 *      Author: babis
 */


#include <iostream>

#include "../hpp/hist.hpp"

void hist( int * x, int * y, int q )
{
	for( int q_i=0; q_i<q; ++q_i )
		y[x[q_i]]++;
}

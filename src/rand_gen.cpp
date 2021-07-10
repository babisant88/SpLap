/*
 * rand_gen.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: babis
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include "../hpp/rand_gen.hpp"

using namespace std;

#define dbg 0

void cumsum( double *, int, double * );
int find( double *, int, double );
// int find( int *, int, int );

/* !!! deprecated !!!
void rand_gen( vector<int>& x, int * y, double * pmf, int num_of_edges, int N )
{
	double * a;

	a = new double[num_of_edges+1];
	a[0] = 0.0;
	cumsum( a, num_of_edges, pmf );

	int * index = new int[N]();

	srand ( time(NULL) );

	int *b = new int[num_of_edges+1];

	for( int n_i=0; n_i<N; ++n_i )
	{
		double r;

		do
			r = ( double )rand()/RAND_MAX;
		while( r == 0 || r == 1.0 );

		double scf = 1.0/r;

		//auto start = chrono::steady_clock::now();
		
		for( int e_i=0; e_i<num_of_edges+1; ++e_i )
		{
			double temp = scf * a[e_i];
			b[e_i] = ( int )( temp > 1.0 );
		}

		index[n_i] = find( b, num_of_edges+1, 1 );
	}

	for( int n_i=0; n_i<N; ++n_i )
		y[n_i] = x[index[n_i]-1];
}
*/



void cumsum( double * a, int num_of_edges, double * pmf )
{
	for( int e_i=1; e_i<num_of_edges+1; ++e_i )
		a[e_i] = a[e_i-1] + pmf[e_i-1];
}



/* !!! deprecated !!!
int find( int *intArray, int size, int value )
{
	int lb = 0;
	int ub = size -1;
	int mid = -1;

	int index = -1;

	while( lb <= ub )
	{
		mid = (ub+lb)/2;

		if( intArray[mid] == value )
		{
			index = mid;
			break;
		}
		else
		{
			// if value is larger
			if( intArray[mid] < value )
			{
				// data is in upper half
				lb = mid + 1;
			}
			// value is smaller
			else
			{
				// data is in lower half
				ub = mid-1;
			}
		}
	}

	while( intArray[index] )
		index = index-1;

	return (index+1);
}
*/



void rand_gen( vector<int>& x, int *y, double * pmf, int num_of_edges, int N )
{
	double * a;

	a = new double[num_of_edges+1];
	a[0] = 0.0;
	cumsum( a, num_of_edges, pmf );

	/* initialize random seed: */
	srand ( time(NULL) );

	for( int n_i=0; n_i<N; ++n_i )
	{
		double r;

		do
			r = ( double )rand()/RAND_MAX;
		while( r == 0 || r == 1.0 );

		double scf = 1.0/r;

		//auto start = chrono::steady_clock::now();

		int index = find(a, num_of_edges+1, 1.0/scf);

//		int e_i;
//		for( e_i=0; e_i<num_of_edges+1; ++e_i )
//		{
//			double temp = scf * a[e_i];
//			if ( temp > 1.0 )
//				break;
//		}

		y[index] += 1;
//		y[e_i-1] += 1;
	}
}



int find( double *arr, int size, double value )
{
	int lb = 0;
	int ub = size -1;
	int mid = -1;

	int index = -1;

	while( lb < ub )
	{
		mid = (ub+lb)/2;

		if( arr[mid] == value )
		{
			index = mid;
			break;
		}
		else
		{
			// if value is larger
			if( arr[mid] < value )
			{
				// data is in upper half
				lb = mid + 1;
			}
			// value is smaller
			else
			{
				// data is in lower half
				ub = mid-1;
			}
		}
	}

	if( lb < ub )
		return index;
	else // lb == ub
	{
		if( arr[lb] > value )
			return (lb-1);
		else
			return lb;
	}
}

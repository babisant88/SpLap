/*
 * rand_gen.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: babis
 */

#include <ctime>       /* time */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include "../hpp/rand_gen.hpp"

using namespace std;

#define dbg 0

void cumsum( double *, int, double * );
int find( int *, int, int );

void rand_gen( vector<int>& x, int * y, double * pmf, int num_of_edges, int N )
{
	double * a;

	a = new double[num_of_edges+1];
	a[0] = 0.0;
	cumsum( a, num_of_edges, pmf );

	int * index = new int[N]();
//	double * c = new double[num_of_edges+1];

	/* initialize random seed: */
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

		index[n_i] = find( b, num_of_edges, 1 );
	
		//auto end = chrono::steady_clock::now();

		//cout << "Elapsed time in microseconds : "
		//<< chrono::duration_cast<chrono::microseconds>(end - start).count()
		//<< " µs" << endl;
	}

	for( int n_i=0; n_i<N; ++n_i )
		y[n_i] = x[index[n_i]-1];
}



void cumsum( double * a, int num_of_edges, double * pmf )
{
	for( int e_i=1; e_i<num_of_edges+1; ++e_i )
		a[e_i] = a[e_i-1] + pmf[e_i-1];
}



int find( int *intArray, int size, int data )
{
	int lowerBound = 0;
	int upperBound = size -1;
	int midPoint = -1;
	int comparisons = 0;
	int index = -1;

	while( lowerBound <= upperBound )
	{

#if dbg
		printf( "Comparison %d\n" , (comparisons +1) );
		printf( "lowerBound : %d, intArray[%d] = %d\n",lowerBound,lowerBound, intArray[lowerBound] );
		printf( "upperBound : %d, intArray[%d] = %d\n",upperBound,upperBound, intArray[upperBound] );
#endif

		comparisons++;

		// compute the mid point
		// midPoint = (lowerBound + upperBound) / 2;
		midPoint = lowerBound + (upperBound - lowerBound) / 2;

		// data found
		if( intArray[midPoint] == data && intArray[midPoint-1] == 0 )
		{
			index = midPoint;
			break;
		}
		else
		{
			// if data is larger
			if( intArray[midPoint] < data )
				// data is in upper half
				lowerBound = midPoint + 1;

			// data is smaller
			else
				// data is in lower half
				upperBound = midPoint -1;
		}
	}

	// printf("Total comparisons made: %d" , comparisons);
	return index;
}

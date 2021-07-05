/*
 * main.cpp
 *
 *  Created on: May 27, 2018
 *      Author: babis
 */

#include <iostream>
#include "../hpp/aux.hpp"
#include "../hpp/graph.hpp"
#include "../hpp/Sparsify.hpp"
#include "../hpp/load_input.hpp"
#include "../hpp/hist.hpp"

#include <omp.h>
#include <chrono>

using namespace std;

#define dbg 0
#define profile 0

int main(int argc, char *argv[])
{
	int n; // matrix size

	n = atoi(argv[2]);

	cout << "n = " << n << endl;

	MatrixXd X(n,n); // Laplacian matrix

#if profile
	auto start = chrono::high_resolution_clock::now();
#endif

	load_input(argv, X);

#if profile
	auto finish = chrono::high_resolution_clock::now();

	chrono::duration<double> elapsed = finish - start;

	cout << "time to load input: " << elapsed.count() << endl;
#endif

	SpMat X_sp;

#if profile
	start = chrono::high_resolution_clock::now();
#endif

	X_sp = Sparsify(X, n, (double)(1/sqrt(n)));

#if profile
	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "time to sparsify laplacian matrix: " << elapsed.count() << endl;
#endif

	cout << "Sparsity ratio of X_sp: " << 1.0-((double)X_sp.nonZeros()/(n*n)) << endl;

	printSpMat( X_sp );
}

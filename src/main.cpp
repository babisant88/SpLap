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

int main(int argc, char *argv[])
{
	int n; // matrix size

	n = atoi(argv[2]);

	cout << "n = " << n << endl;

	MatrixXd X(n, n); // Laplacian matrix

	auto start = chrono::high_resolution_clock::now();
	load_input(argv, X);
	auto finish = chrono::high_resolution_clock::now();

	chrono::duration<double> elapsed = finish - start;

	cout << "time to load input: " << elapsed.count() << endl;

	SpMat X_sp;

	start = chrono::high_resolution_clock::now();
	X_sp = Sparsify(X, n, (double)(1/sqrt(n)));

	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "time to sparsify laplacian matrix: " << elapsed.count() << endl;

	cout << "Sparsity ratio of X_sp: " << 1.0-((double)X_sp.nonZeros()/(n*n)) << endl;

	//printSpMatToFile(X_sp, "X_sp.dat");
}

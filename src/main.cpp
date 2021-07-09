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
#define reg_test 1

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

	X_sp = Sparsify_top(X, n);

#if profile
	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "time to sparsify laplacian matrix: " << elapsed.count() << endl;
#endif

	cout << "Sparsity ratio of X_sp: " << 1.0-((double)X_sp.nonZeros()/(n*n)) << endl;

#if reg_test

	VectorXd b = VectorXd::Random(n);
	VectorXd golden_sol, sol, diff;

	LLT<MatrixXd> solver_d;
	solver_d.compute(X);
	golden_sol = solver_d.solve(b);

	/* for a fair (due to different solvers (dense/sparse)) comparison,
	 * convert sparse matrix X_sp to Dense */
	MatrixXd X_d = MatrixXd(X_sp);

	LLT<MatrixXd> solver_sp;
	solver_sp.compute(X_d);
	sol = solver_sp.solve(b);

	diff = sol - golden_sol;

	cout << "relative error = " << (diff.norm())/(golden_sol.norm()) << endl;

#endif
}

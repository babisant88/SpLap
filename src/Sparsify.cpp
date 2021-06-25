/*
 * Sparsify.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: babis
 */

#include <cmath>
#include <iostream>

//#include "../hpp/aux.hpp"
#include "../hpp/graph.hpp"
#include "../hpp/pinv.hpp"
#include "../hpp/Sparsify.hpp"
#include "../hpp/rand_gen.hpp"
#include "../hpp/hist.hpp"

#include <chrono>

#define dbg 0
#define rand_smp 0

using namespace std;

void find(SpMat&, double*, int*);
bool check_Cval(double, double, int);

/* Precondition: the input matrix M has to be Symmetric */
SpMat Sparsify(MatrixXd& M, int n, double epsilon)
{
	SpMat M_sp;

	double *pe, *Reff;

	MatrixXd L, A;
	MatrixXd pinv_L; // the pseudo-inverse of the Laplacian matrix
	VectorXd D, res_to_grd;

	graph * A_grph;

	/* MATLAB: A = abs(diag(diag(M)) - M); */
	A = M.diagonal().asDiagonal(); // A is the Adjacency matrix of A_grph
	A = A - M;
	A = A.cwiseAbs();

#if dbg
	cout << "A = " << endl;
	printMatrixXd(A);
#endif

	/* MATLAB: D = sum(A); */
	D = A.rowwise().sum(); // D is the Degree matrix of A_grph

#if dbg
	cout << "D = " << endl;
	printVectorXd(D);
#endif

	/* MATLAB: res_to_grd = diag(M) - D'; */
	res_to_grd = M.diagonal();
	res_to_grd = res_to_grd - D;

#if dbg
	cout << "res_to_grd = " << endl;
	printVectorXd(res_to_grd);
#endif

	/* MATLAB: L = diag(D) - A; */
	L = D.asDiagonal(); // L is the Laplacian matrix of A_grph
	L = L - A;

#if dbg
	cout << "L = " << endl;
	printMatrixXd(L);
#endif

	/* MATLAB: A_grph = graph(A); */
	A_grph = new graph(A);

	auto start = chrono::high_resolution_clock::now();

	/* Compute L^(+) */
	pinv_L = pinv<MatrixXd>(L);

#if dbg
	cout << "L^(+) = " << endl;
	printMatrixXd(pinv_of_L);
#endif

	Reff = new double[A_grph->num_of_edges];

	SpMat B = A_grph->get_incidence_matrix();

	/* Obtain values & innerIndices of B */
	double * B_x = new double[2*A_grph->num_of_edges];
	int * B_i = new int[2*A_grph->num_of_edges];

	find(B, B_x, B_i);

	for(int e_i=0; e_i<A_grph->num_of_edges; ++e_i)
	{
		double tmp;

		tmp = B_x[2*e_i] * ( B_x[2*e_i]*pinv_L(B_i[2*e_i], B_i[2*e_i]) + B_x[2*e_i+1]*pinv_L(B_i[2*e_i+1], B_i[2*e_i]) );
		Reff[e_i] = tmp + B_x[2*e_i+1] * ( B_x[2*e_i]*pinv_L(B_i[2*e_i], B_i[2*e_i+1]) + B_x[2*e_i+1]*pinv_L(B_i[2*e_i+1], B_i[2*e_i+1]) );
	}

	delete[] B_x;
	delete[] B_i;

	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;

	cout << "time to compute Reff: " << elapsed.count() << endl;

#if dbg
	cout << "Reff = " << endl;
	for(int row_i=0; row_i<A_grph->num_of_edges; row_i++)
		cout << Reff[row_i] << '\t';

	cout << endl << endl;
#endif

	/* define a PMF over the set of edges as follows: */
	/* for each e: p(e) = (w(e)*Reff(e))/sum[for each e: w(e)*Reff(e)] */

	double sum = 0.0;
	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
		sum += ( A_grph->we[e_i] * Reff[e_i] );

	pe = new double[A_grph->num_of_edges];

	for( int pe_i=0; pe_i<A_grph->num_of_edges; ++pe_i )
		pe[pe_i] = ( A_grph->we[pe_i] * Reff[pe_i] )/( sum );

#if dbg
	sum = 0.0;
	for( int pe_i=0; pe_i<A_grph->num_of_edges; ++pe_i )
		sum += pe[pe_i];

	cout << "n = " << n << endl;
	cout << "sum[pe] = " << sum << endl;
#endif

	delete[] Reff;

	const double C = 0.1;

	bool valid_Cval = check_Cval(epsilon, C, n);

	cout << "valid_Cval = " << valid_Cval << endl;

	if ( !valid_Cval )
	{
		cout << "C = " << C << " is an unacceptable value" << endl;
		exit(-1);
	}

	/* number of Monte Carlo trials */
	double q = ( 9*(C*C)*n*log(n) )/( epsilon*epsilon );

#if dbg
	cout << "q = " << (int)ceil(q) << endl;
#endif

#if rand_smp
	int * H_edges_i;

	H_edges_i = new int[(int)ceil(q)];

	rand_gen( A_grph->ide, H_edges_i, pe, A_grph->num_of_edges, (int)ceil(q) );

	int * hist_occ = new int[A_grph->num_of_edges]();
	hist( H_edges_i, hist_occ, (int)ceil(q) );

	delete[] H_edges_i;
#else
	/* how many times an edge has been chosen? */
	int * hist_occ = new int[A_grph->num_of_edges];

	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
		hist_occ[e_i] = ( int )round( pe[e_i] * q );
#endif

#if dbg
	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
		cout << "hist_occ[" << e_i << "] = " << hist_occ[e_i] << endl;
#endif

	double * H_edges_wght = new double[A_grph->num_of_edges];

	/* set the weights of edges to be inserted in H_grph */
	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
		H_edges_wght[e_i] = ( A_grph->we[e_i]*hist_occ[e_i] )/( ceil(q)*pe[e_i] );

	delete[] pe;
	delete[] hist_occ;

	graph * H_grph = new graph(n);

	start = chrono::high_resolution_clock::now();

	/* add to H_grph the selected edges (those edges that have non-zero weight)*/
	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
		if( H_edges_wght[e_i] != 0 )
			H_grph->add_edge( A_grph->vertex0e[e_i], A_grph->vertex1e[e_i], H_edges_wght[e_i] );

	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "time to build H_grph: " << elapsed.count() << endl;

	delete A_grph;
	delete[] H_edges_wght;

	SpMat H_lpc;
	H_lpc = H_grph->get_laplacian_matrix();

	/* MATLAB: M_sp = sparse(diag(res_to_grd)) + H_lpc; */
	M_sp = res_to_grd.asDiagonal();
	M_sp += H_lpc;

	delete H_grph;

	return M_sp;
}

void find(SpMat& B, double * B_x, int * B_i)
{
	int nz_i = 0;

	for(int k=0; k<B.outerSize(); ++k)
		for(SpMat::InnerIterator it(B, k); it; ++it)
		{
			B_x[nz_i] = it.value();
			B_i[nz_i] = it.index(); // inner index
			nz_i++;
		}
}

bool check_Cval(double epsilon, double C, int n)
{
	double S;
	bool valid;

	int M = n-1;

	S = C * sqrt( (epsilon*epsilon) * (( log( ( 9*(C*C)*n*log(n) )/(epsilon*epsilon) ) * ( M ) )/( 9*(C*C)*n*log(n) )) );

#if dbg
	cout << "epsilon = " << epsilon << endl;
	cout << "S = " << S << endl;
#endif

	if ( S <= epsilon/2 )
	    valid = true;
	else
	    valid = false;

	return valid;
}

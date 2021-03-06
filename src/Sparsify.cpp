/*
 * Sparsify.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: babis
 */

#include <cmath>
#include <iostream>
#include <chrono>

#include "../hpp/aux.hpp"
#include "../hpp/graph.hpp"
#include "../hpp/pinv.hpp"
#include "../hpp/Sparsify.hpp"
#include "../hpp/rand_gen.hpp"
#include "../hpp/hist.hpp"
#include "../hpp/ApproxReff.hpp"

#define dbg 0 // enable debugging info
#define rand_smp 0 // enable random sampling
#define profile 0 // endalbe profiling info
#define Reffx 0 // enable the approximation of the effective resistances

using namespace std;

bool check_Cval(double, double, int);



SpMat Sparsify_top(MatrixXd& M, int n)
{
	SpMat M_sp;

	MatrixXd L, A;

	VectorXd D, res_to_grd;

	graph * A_grph;

	/* A = abs(diag(diag(M)) - M); */
	A = M.diagonal().asDiagonal(); // A is the Adjacency matrix of A_grph
	A = A - M;
	A = A.cwiseAbs();

#if dbg
	cout << "A = " << endl;
	printMatrixXd(A);
#endif

	/* D = sum(A); */
	D = A.rowwise().sum(); // D is the Degree matrix of A_grph

#if dbg
	cout << "D = " << endl;
	printVectorXd(D);
#endif

	/* res_to_grd = diag(M) - D'; */
	res_to_grd = M.diagonal();
	res_to_grd = res_to_grd - D;

#if dbg
	cout << "res_to_grd = " << endl;
	printVectorXd(res_to_grd);
#endif

	/* L = diag(D) - A; */
	L = D.asDiagonal(); // L is the Laplacian matrix of A_grph
	L = L - A;

#if dbg
	cout << "L = " << endl;
	printMatrixXd(L);
#endif

	/* A_grph = graph(A); */
	A_grph = new graph(L);

	/* Find the Connected Components (CCs) of graph A_grph */
	vector<vector<int>> CCs = A_grph->ConnComp();

#if dbg
	for( uint i=0; i<CCs.size(); ++i )
	{
		for( uint j=0; j<CCs[i].size(); ++j )
			cout << CCs[i][j] << " ";
		cout << endl;
	}
#endif

	/* initialize the resulting graph (hopefully sparse) */
	graph * H_grph = new graph(n);

	/* loop over the CCs */
	for( uint CC_i=0; CC_i<CCs.size(); ++CC_i )
	{
		uint CCi_n = CCs[CC_i].size(); // number of nodes in CCi

		/* Extract the Laplacian matrix of the i-th Connected Component (CCi) */
		MatrixXd CCi_L_tmp(CCi_n,n);

		for( uint v_i=0; v_i<CCi_n; ++v_i )
			CCi_L_tmp.row(v_i) = L.row( CCs[CC_i][v_i] );

		MatrixXd CCi_L(CCi_n, CCi_n);

		for( uint v_i=0; v_i<CCi_n; ++v_i )
			CCi_L.col(v_i) = CCi_L_tmp.col( CCs[CC_i][v_i] );

#if dbg
		printMatrixXd(CCi_L);
#endif
		/* mapping between the edges of the initial graph and
		 * the edges of the i-th Connected Component */
		vector<int> edge_map;

		uint CCi_e= 0;

		for( uint v_i=0; v_i<CCi_n; ++v_i )
			for( list<int>::iterator iter = A_grph->adj[CCs[CC_i][v_i]].begin(); iter != A_grph->adj[CCs[CC_i][v_i]].end(); ++iter )
			{
				edge_map.push_back(*iter);
				CCi_e++;
			}
#if dbg
		for( uint i=0; i<CCi_e; ++i )
			cout << i << " -> " << edge_map[i] << endl;
#endif

#if profile
		auto start = chrono::high_resolution_clock::now();
#endif

		graph *CCi_grph = NULL;

		double *H_edges_w = Sparsify(CCi_L, CCi_n, CCi_grph, false);

#if dbg
		for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
			cout << "H_edges_w[" << e_i << "] = " << H_edges_w[e_i] << endl;
#endif

		/* add to H_grph the selected edges (those edges that have non-zero weight) */
		for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
			if( H_edges_w[e_i] != 0 )
				H_grph->add_edge( A_grph->vertex0e[edge_map[e_i]], A_grph->vertex1e[edge_map[e_i]], H_edges_w[e_i] );

#if profile
		auto finish = chrono::high_resolution_clock::now();
		chrono::duration<double> elapsed = finish - start;
		cout << "time to build H_grph: " << elapsed.count() << endl;
#endif

		delete CCi_grph;
		delete[] H_edges_w;
	}

	SpMat H_lpc;
	H_lpc = H_grph->get_laplacian_matrix_sp();

	/* M_sp = sparse(diag(res_to_grd)) + H_lpc; */
	M_sp = res_to_grd.asDiagonal();
	M_sp += H_lpc;

	delete A_grph;
	delete H_grph;

	return M_sp;
}



double* Sparsify(MatrixXd& CCi_L, int CCi_n, graph*& CCi_grph, bool MLST_en)
{
	double *pe;

	/* error/sparsity trade-off parameter */
	/* epsilon has to be between 1/sqrt(CCi_n) and 1 */
	/* larger values of epsilon lead to more sparse but less accurate model */
	double epsilon = 0.4; //1/sqrt((double)CCi_n);

	MatrixXd pinv_L; // the pseudo-inverse of the Laplacian matrix

	CCi_grph = new graph(CCi_L);

#if profile
	auto start = chrono::high_resolution_clock::now();
#endif

#if !Reffx
	/* Compute L^(+) */
	pinv_L = pinv<MatrixXd>(CCi_L);
#endif

#if dbg
	cout << "L^(+) = " << endl;
	printMatrixXd(pinv_L);
#endif

#if Reffx

	double *Reff = appxReff_mem(CCi_grph, 0.1);

#else
	double *Reff = new double[CCi_grph->num_of_edges];

	SpMat B = CCi_grph->get_incidence_matrix();

	/* Obtain values & innerIndices of B */
	double * B_x = new double[2*CCi_grph->num_of_edges];
	int * B_i = new int[2*CCi_grph->num_of_edges];

	find(B, B_x, B_i);

	for(int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i)
	{
		double tmp;

		tmp = B_x[2*e_i] * ( B_x[2*e_i]*pinv_L(B_i[2*e_i], B_i[2*e_i]) + B_x[2*e_i+1]*pinv_L(B_i[2*e_i+1], B_i[2*e_i]) );
		Reff[e_i] = tmp + B_x[2*e_i+1] * ( B_x[2*e_i]*pinv_L(B_i[2*e_i], B_i[2*e_i+1]) + B_x[2*e_i+1]*pinv_L(B_i[2*e_i+1], B_i[2*e_i+1]) );
	}

#if profile
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;

	cout << "time to compute Reff: " << elapsed.count() << endl;
#endif

	delete[] B_x;
	delete[] B_i;

#endif

#if dbg
	cout << "Reff = " << endl;
	for(int row_i=0; row_i<CCi_grph->num_of_edges; row_i++)
		cout << Reff[row_i] << '\t';

	cout << endl << endl;

#endif

	/* define a PMF over the set of edges as follows: */
	/* for each e: p(e) = (w(e)*Reff(e))/sum[for each e: w(e)*Reff(e)] */
	double sum = 0.0;
	for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
		sum += ( CCi_grph->we[e_i] * Reff[e_i] );

	pe = new double[CCi_grph->num_of_edges];

	for( int pe_i=0; pe_i<CCi_grph->num_of_edges; ++pe_i )
		pe[pe_i] = ( CCi_grph->we[pe_i] * Reff[pe_i] )/( sum );

#if dbg
	for( int pe_i=0; pe_i<CCi_grph->num_of_edges; ++pe_i )
		cout << "pe[" << pe_i << "] = " << pe[pe_i] << endl;

	sum = 0.0;
	for( int pe_i=0; pe_i<CCi_grph->num_of_edges; ++pe_i )
		sum += pe[pe_i];

	cout << "CCi_n = " << CCi_n << endl;
	cout << "sum[pe] = " << sum << endl;
#endif

	double C = 1.0;
	/* loop until you find a valid value for C */
	while( !check_Cval(epsilon, C, CCi_n) )
		C /= 2;

#if dbg
	cout << "C = " << C << endl;
#endif

	/* number of Monte Carlo trials */
	double q = ( 9*(C*C)*CCi_n*log(CCi_n) )/( epsilon*epsilon );

	vector<int> MLSTedges;

	if ( MLST_en )
	{
		/* Compute the "Maximum Likelihood" Spanning Tree (MLST) to  */
		/* avoid a disjoint graph as a result after sparsification */
		graph *CCi_likelihood_grph = new graph(CCi_L);

		for( int e_i = 0; e_i<CCi_likelihood_grph->num_of_edges; ++e_i )
			CCi_likelihood_grph->we[e_i] = 1.0/pe[e_i];

		MLSTedges = CCi_likelihood_grph->MST(); // MLSTedges are the edges that belong to MLST

		delete CCi_likelihood_grph;

#if dbg
	for( std::vector<int>::iterator e_it = MLSTedges.begin(); e_it != MLSTedges.end(); e_it++ )
		cout << (*e_it) << " ";

	cout << endl;
#endif

		/* Subtrack from q the number of edges belonging to the MLST */
		q -= (MLSTedges.size()-1);
	}

#if dbg
	cout << "q = " << (int)ceil(q) << endl;
#endif

	delete[] Reff;

	double *H_edges_w = new double[CCi_grph->num_of_edges]();

#if rand_smp

/* !!! deprecated !!!
	int * H_edges_i;

	H_edges_i = new int[(int)ceil(q)];

	rand_gen( CCi_grph->ide, H_edges_i, pe, CCi_grph->num_of_edges, (int)ceil(q) );

#if dbg
	for( int t_i=0; t_i<(int)ceil(q); ++t_i )
		cout << "H_edges_i[" << t_i << "] = " << H_edges_i[t_i] << endl;
#endif

	int * hist_occ = new int[CCi_grph->num_of_edges]();
	hist( H_edges_i, hist_occ, (int)ceil(q) );

	delete[] H_edges_i;
*/

	int * hist_occ = new int[CCi_grph->num_of_edges]();
	rand_gen( CCi_grph->ide, hist_occ, pe, CCi_grph->num_of_edges, (int)ceil(q) );

#if dbg
	for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
		cout << "hist_occ[" << e_i << "] = " << hist_occ[e_i] << endl;
#endif

	/* add the MLST edges for sure (recall: to avoid disjoint graph)*/
	for( uint e_i=0; e_i < MLSTedges.size(); ++e_i )
		hist_occ[MLSTedges[e_i]] += 1;

	/* set the weights of edges to be inserted in H_grph */
	for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
		H_edges_w[e_i] = ( CCi_grph->we[e_i]*hist_occ[e_i] )/( ceil(q)*pe[e_i] );

	delete[] hist_occ;
#else

	/* Pick the edges to be inserted in H_grph */
	for( int e_i=0; e_i<CCi_grph->num_of_edges; ++e_i )
		/* I expect the condition below to return 1 in case it is true and 0 otherwise */
		H_edges_w[e_i] = CCi_grph->we[e_i] * (( (int)round( pe[e_i]*q ) ) > 0);

	/* add the MLST edges for sure (recall: to avoid disjoint graph)*/
	for( uint e_i=0; e_i < MLSTedges.size(); ++e_i )
		H_edges_w[MLSTedges[e_i]] = CCi_grph->we[e_i];
#endif

	delete[] pe;

	return H_edges_w;
}



bool check_Cval(double epsilon, double C, int n)
{
	double S;
	bool valid;

	int M = n-1;

	double _log_ = log( ( 9*(C*C)*n*log(n) )/(epsilon*epsilon) );

	S = C * sqrt( (epsilon*epsilon) * (( _log_  * ( M ) )/( 9*(C*C)*n*log(n) )) );

#if dbg
	cout << "epsilon = " << epsilon << endl;
	cout << "S = " << S << endl;
#endif

	if ( (_log_ >= 0) && ( S <= epsilon/2 ) )
	    valid = true;
	else
	    valid = false;

	return valid;
}

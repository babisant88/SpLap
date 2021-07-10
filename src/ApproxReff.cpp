/*
 * ApproxReff.cpp
 *
 *  Created on: Jul 3, 2021
 *      Author: babis
 */

#include "../hpp/ApproxReff.hpp"
#include <iostream>
#include <cmath>
#include "../hpp/iter_solver.hpp"

#define dbg 0


double *appxReff_mem(graph *A_grph, double delta)
{
	double *Reff = new double[A_grph->num_of_edges]();

	uint k = (uint)ceil(log(A_grph->num_of_vertices)/(delta*delta));

	SpMat B = A_grph->get_incidence_matrix();

#if dbg
	cout << "B = " << endl;
	printSpMat(B);
#endif

	/* Obtain values & innerIndices of B */
	double *B_x = new double[2*A_grph->num_of_edges];
	int *B_i = new int[2*A_grph->num_of_edges];

	find(B, B_x, B_i);

	for( int e_i=0; e_i < A_grph->num_of_edges; ++e_i )
	{
		B_x[2*e_i] *= sqrt(fabs(A_grph->we[e_i]));
		B_x[2*e_i+1] *= sqrt(fabs(A_grph->we[e_i]));
	}

	SpMat Ag(A_grph->num_of_edges,A_grph->num_of_vertices);

	std::vector<T> nz;

	/* Ag = W^(1/2) * B^(T) */

	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
	{
		nz.push_back(T(e_i, B_i[2*e_i], B_x[2*e_i]));
		nz.push_back(T(e_i, B_i[2*e_i+1], B_x[2*e_i+1]));
	}

	Ag.setFromTriplets(nz.begin(),nz.end());

#if dbg
	cout << "Ag = " << endl;
	printSpMat(Ag);
#endif

	MatrixXd Lg = A_grph->get_laplacian_matrix_d();

#if dbg
	cout << "Lg = " << endl;
	printMatrixXd(Lg);
#endif

	LLT<MatrixXd> solver;
	solver.compute( Lg ); // Lg is expected to be dense...that's why we want its sparsification

/*	In case Lg matrix is sparse, I suggest you to try the cg solver
 *  (with a suitable prconditioner e.g.. CMG from Koutis, Miller & Peng */

//	ConjugateGradient<MatrixXd, Lower|Upper> cg;
//	cg.compute(Lg);

	for (uint row_i=0; row_i<k; ++row_i)
	{
		VectorXd Z;

		/* A vector with random values between -1 and +1*/
		VectorXd Q = VectorXd::Random(A_grph->num_of_edges);

		/* round Qi entries to the nearest nozero integer value (+/-)1 */
		for( int q_i=0; q_i<Q.size(); q_i++ )
		{
				if( Q(q_i) < 0 )
					Q(q_i) = -1;
				else
					Q(q_i) = +1;
		}

#if dbg
		cout << "Q = " << endl;
		printVectorXd(Q);
#endif

		Q = 1/sqrt((double)k) * Q;

		VectorXd Y = Q.transpose()*Ag;

#if dbg
		cout << "Y = " << endl;
		printVectorXd(Y);
#endif

		Z = solver.solve(Y);
//		Z = cg.solve(b);

		for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
			Reff[e_i] += pow( Z(A_grph->vertex1e[e_i]) - Z(A_grph->vertex0e[e_i]), 2.0 );
	}

	return Reff;
}



double *appxReff(graph *A_grph, double delta)
{
	double *Reff = new double[A_grph->num_of_edges];

	uint k = (uint)ceil(log(A_grph->num_of_vertices)/(0.1*0.1));

	/* A matrix with random values between -1 and +1*/
	MatrixXd Q = MatrixXd::Random(k,A_grph->num_of_edges);

	/* round Q entries to the nearest nozero integer value (+/-)1 */
	for( int row_i=0; row_i<Q.rows(); row_i++ )
		for( int col_i=0; col_i<Q.cols(); col_i++ )
		{
			if( Q(row_i,col_i) < 0 )
				Q(row_i, col_i) = -1;
			else
				Q(row_i, col_i) = +1;
		}

#if dbg
	cout << "Q = " << endl;
	printMatrixXd(Q);
#endif

	Q = 1/sqrt((double)k) * Q;

#if dbg
	cout << "Q = " << endl;
	printMatrixXd(Q);
#endif

	SpMat B = A_grph->get_incidence_matrix();

#if dbg
	cout << "B = " << endl;
	printSpMat(B);
#endif

	/* Obtain values & innerIndices of B */
	double *B_x = new double[2*A_grph->num_of_edges];
	int *B_i = new int[2*A_grph->num_of_edges];

	find(B, B_x, B_i);

	/* Ag = W^(1/2) * B^(T) */

	for( int e_i=0; e_i < A_grph->num_of_edges; ++e_i )
	{
		B_x[2*e_i] *= sqrt(fabs(A_grph->we[e_i]));
		B_x[2*e_i+1] *= sqrt(fabs(A_grph->we[e_i]));
	}

	SpMat Ag(A_grph->num_of_edges,A_grph->num_of_vertices);

	std::vector<T> nz;

	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
	{
		nz.push_back(T(e_i, B_i[2*e_i], B_x[2*e_i]));
		nz.push_back(T(e_i, B_i[2*e_i+1], B_x[2*e_i+1]));
	}

	Ag.setFromTriplets(nz.begin(),nz.end());

#if dbg
		cout << "Ag = " << endl;
		printSpMat(Ag);
#endif

	MatrixXd Y(k,A_grph->num_of_vertices);

	Y = Q*Ag;

#if dbg
		cout << "Y = " << endl;
		printMatrixXd(Y);
#endif

	MatrixXd Lg = A_grph->get_laplacian_matrix_d();

#if dbg
		cout << "Lg = " << endl;
		printMatrixXd(Lg);
#endif

	LLT<MatrixXd> solver;

	solver.compute(Lg); // Lg is expected to be dence...that's why we want its sparsification

	/*	In case Lg matrix is sparse, I suggest you to try the cg solver
	 *  (with a suitable prconditioner e.g.. CMG from Koutis, Miller & Peng */

//	ConjugateGradient<MatrixXd, Lower|Upper> cg;
//	cg.compute(Lg);

	MatrixXd Z(k,A_grph->num_of_vertices);

	//VectorXd x = VectorXd::Zero(Lg.rows());

	/* Z(i,:) = LapSolve(L, Y(i,:)) */
	for (int row_i=0; row_i<Y.rows(); ++row_i)
	{
		VectorXd b = Y.row(row_i);
//		Z.row(row_i) = cg.solve(b);

		Z.row(row_i) = solver.solve(b);
	}

	/* Reff(e) = ||Z * ( e(i) - e(j) )||^2 */
	for( int e_i=0; e_i<A_grph->num_of_edges; ++e_i )
	{
		VectorXd re = Z.col( A_grph->vertex1e[e_i] ) - Z.col( A_grph->vertex0e[e_i] );
		Reff[e_i] = pow(re.norm(), 2.0);
	}

	return Reff;
}

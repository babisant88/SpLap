/*
 * graph.cpp
 *
 *  Created on: May 29, 2018
 *      Author: babis
 */

#include <iostream>

#include "../hpp/graph.hpp"

// initialize the graph with its adjacency matrix
graph::graph(MatrixXd& A)
{
	int n = A.rows();

	num_of_vertices = n;

	int ei = 0;

	for(int row_i=0; row_i<n-1; row_i++)
	{
		for(int col_i=row_i+1; col_i<n; col_i++)
		{
			if ( A(row_i, col_i) != 0 )
			{
				ide.push_back( ei++ );
				vertex0e.push_back( row_i );
				vertex1e.push_back( col_i );
				we.push_back( A(row_i, col_i) );
			}
		}
	}

	num_of_edges = ei;
}

graph::graph( int n )
{
	num_of_vertices = n;
	num_of_edges = 0;
}

graph::~graph()
{
}

void graph::add_edge(int i, int j, double w)
{
	ide.push_back( num_of_edges++ );
	vertex0e.push_back( i );
	vertex1e.push_back( j );
	we.push_back( w );
}

SpMat graph::get_laplacian_matrix()
{
	SpMat L(this->num_of_vertices, this->num_of_vertices);

	std::vector<T> nz;

	if ( num_of_edges != 0 )
	{
		for( int ei=0; ei<num_of_edges; ++ei )
		{
			nz.push_back(T(vertex0e[ei], vertex1e[ei], -we[ei]));
			nz.push_back(T(vertex1e[ei], vertex0e[ei], -we[ei]));
		}

		double * row_sums = new double[num_of_vertices]();

		for( int e_i=0; e_i<2*num_of_edges; ++e_i )
			row_sums[nz[e_i].row()] += abs( nz[e_i].value() );

		for( int r_i=0; r_i<num_of_vertices; ++r_i )
			nz.push_back( T(r_i, r_i, row_sums[r_i]) );

		L.setFromTriplets(nz.begin(), nz.end());
	}
	else
		cout << "Be careful! L matrix is empty...graph contains 0 edges" << endl;

	return L;
}

SpMat graph::get_incidence_matrix()
{
	SpMat B(this->num_of_vertices, this->num_of_edges);

	std::vector<T> nz;

	if ( num_of_edges != 0 )
	{
		for( int ei=0; ei<num_of_edges; ++ei )
		{
			nz.push_back(T(vertex0e[ei], ei, 1.0));
			nz.push_back(T(vertex1e[ei], ei, -1.0));
		}

		B.setFromTriplets(nz.begin(), nz.end());
	}
	else
		cout << "Be careful! B matrix is empty...graph contains 0 edges" << endl;

	return B;
}

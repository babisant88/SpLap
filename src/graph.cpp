/*
 * graph.cpp
 *
 *  Created on: May 29, 2018
 *      Author: babis
 */

#include <iostream>
#include <vector>
#include <queue>

#include "../hpp/graph.hpp"

// initialize the graph with its Laplacian matrix
graph::graph(MatrixXd& L)
{
	int n = L.rows();

	num_of_vertices = n;

	adj = new list<int>[n];

	int ei = 0;

	for(int row_i=0; row_i<n-1; row_i++)
	{
		for(int col_i=row_i+1; col_i<n; col_i++)
		{
			if ( L(row_i, col_i) != 0 )
			{
				adj[row_i].push_back( ei );

				ide.push_back( ei++ );
				vertex0e.push_back( row_i );
				vertex1e.push_back( col_i );
				we.push_back( -L(row_i, col_i) );
			}
		}
	}

	num_of_edges = ei;
}



graph::graph( int n )
{
	num_of_vertices = n;
	num_of_edges = 0;

	adj = new list<int>[n];
}



graph::~graph()
{
	delete[] adj;
}



void graph::add_edge(int i, int j, double w)
{
	adj[i].push_back( num_of_edges );

	ide.push_back( num_of_edges++ );
	vertex0e.push_back( i );
	vertex1e.push_back( j );
	we.push_back( w );
}



SpMat graph::get_laplacian_matrix_sp()
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

		double *row_sums = new double[num_of_vertices]();

		for( int e_i=0; e_i<2*num_of_edges; ++e_i )
			row_sums[nz[e_i].row()] += abs( nz[e_i].value() );

		for( int r_i=0; r_i<num_of_vertices; ++r_i )
			nz.push_back( T(r_i, r_i, row_sums[r_i]) );

		L.setFromTriplets(nz.begin(), nz.end());
		delete[] row_sums;
	}
	else
		cout << "Be careful! L matrix is empty...graph contains 0 edges" << endl;

	return L;
}



MatrixXd graph::get_laplacian_matrix_d()
{
	MatrixXd L = MatrixXd::Zero(num_of_vertices, num_of_vertices);

	if ( num_of_edges != 0 )
	{
		for( int e_i=0; e_i<num_of_edges; ++e_i )
		{
			L(vertex0e[e_i], vertex1e[e_i]) = -we[e_i];
			L(vertex1e[e_i], vertex0e[e_i]) = -we[e_i];
		}

		for( int v_i=0; v_i<num_of_vertices; ++v_i )
		{
			double row_sum = 0;
			for( int v_j=0; v_j<num_of_vertices; ++v_j )
				row_sum += (-L(v_i,v_j));

			L(v_i,v_i) = row_sum;
		}
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



vector<vector<int>> graph::ConnComp()
{
	vector<vector<int>> CCs; // a vector with the CCs

	// initialize vertices as not visited
	vector<bool> visited(num_of_vertices, false);

	queue<int> Q;

	for( int v=0; v < num_of_vertices; v++ )
	{
		vector<int> CCi; // the node indices that belong to the i-th CC

		if ( !visited[v] )
		{
			/* BF traversal for finding the CCs */

			Q.push(v);
			visited[v] = true;

			CCi.push_back(v);

			while( !Q.empty() )
			{
				int n = Q.front();
				Q.pop();

				list<int>::iterator iter;

				for( iter = adj[n].begin(); iter != adj[n].end(); ++iter )
				{
					if( !visited[vertex1e[*iter]] )
					{
						Q.push(vertex1e[*iter]);
						visited[vertex1e[*iter]] = true;

						CCi.push_back(vertex1e[*iter]);
					}
				}
			}

			CCs.push_back( CCi );
		}
	}

	return CCs;
}



/* Kruskal's MST algorithm */
vector<int> graph::MST()
{
	vector<int> MLSTedges;

	vector<int> parent(num_of_vertices);

	vector<int> cluster_size(num_of_vertices, 0);

	for( int v_i=0; v_i<num_of_vertices; ++v_i )
		parent[v_i] = v_i;

	/* Sort all edges into ascending
	 * order by weight w */
	pair<double, int> e_w_pairs[num_of_edges];

	for( int e_i = 0; e_i < num_of_edges; e_i++ )
    {
		e_w_pairs[e_i].first = we[e_i];
		e_w_pairs[e_i].second = e_i;
    }

	sort(e_w_pairs, e_w_pairs + num_of_edges);

	/* for each (u,v) taken from the sorted
	 * list of edges */
	for( int e_i = 0; e_i < num_of_edges; e_i++ )
	{
		int e_j = e_w_pairs[e_i].second;

		int v_i = vertex0e[e_j];
		int v_j = vertex1e[e_j];

		while( parent[v_i] != v_i ) // FIND_SET(v_i): See whether v_i and v_j fall into the same cluster of nodes
			v_i = parent[v_i];

		while( parent[v_j] != v_j ) // FIND_SET(v_j)
			v_j = parent[v_j];

		if( v_i != v_j )
		{
			MLSTedges.push_back(e_j);

			if( cluster_size[v_i] >= cluster_size[v_j] )
			{
				parent[v_j] = v_i; // union(v_i,v_j): Unify the cluster of node v_i belong to with the cluster of nodes v_j belogs to
				cluster_size[v_i] += ( cluster_size[v_j] + 1);
			}
			else
			{
				parent[v_i] = v_j; // union(v_j,v_i)
				cluster_size[v_j] += ( cluster_size[v_i] + 1);
			}
		}
	}

	return MLSTedges;
}

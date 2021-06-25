/*
 * graph.hpp
 *
 *  Created on: May 29, 2018
 *      Author: babis
 */

#ifndef MYHPPFILES_GRAPH_HPP_
#define MYHPPFILES_GRAPH_HPP_

#include "types.hpp"

#include <vector>

using namespace std;

class graph
{
	public:
			graph(MatrixXd &);
			graph( int );
			~graph();

			int num_of_vertices;
			int num_of_edges;

			/* set of edges */
			vector<int> ide;
			vector<int> vertex0e;
			vector<int> vertex1e;
			vector<double> we;

			void add_edge(int, int, double);

			SpMat get_incidence_matrix();
			SpMat get_laplacian_matrix();
};

#endif /* MYHPPFILES_GRAPH_HPP_ */

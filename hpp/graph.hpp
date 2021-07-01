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
#include <list>

using namespace std;

class graph
{
	public:
			graph(MatrixXd &);
			graph( int );
			~graph();

			int num_of_vertices;
			int num_of_edges;

			/* adj list */
			list<int> * adj;

			/* set of edges */
			vector<int> ide;
			vector<int> vertex0e;
			vector<int> vertex1e;
			vector<double> we;

			void add_edge(int, int, double);

			SpMat get_incidence_matrix();
			SpMat get_laplacian_matrix();

			vector<vector<int>> ConnComp();
};

#endif /* MYHPPFILES_GRAPH_HPP_ */

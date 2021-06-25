/*
 * load_input.cpp
 *
 *  Created on: Jun 12, 2018
 *      Author: babis
 */

#include <iostream>
#include <fstream>
#include <cstring>

//#include "../hpp/aux.hpp"
#include "../hpp/load_input.hpp"

#define dbg 0

using namespace std;


void load_input(char* argv[], MatrixXd& X)
{
	load_dense_matrix(argv[1], X);
}

void load_dense_matrix(char* filename, MatrixXd& X)
{
	string line;
	ifstream myfile ( filename );

	int row_i, col_i;

	row_i = 0;

	if ( myfile.is_open() )
	{
		while( getline( myfile, line ) )
		{
			char * pch;
			pch = strtok( (char *)line.c_str(), " " );

			col_i = 0;

			while ( pch != NULL )
			{
				X(row_i, col_i) = atof( pch );
				col_i++;
				pch = strtok ( NULL, " " );
			}

			row_i++;
		}

		myfile.close();
	}
	else
		cout << "Unable to open file" << filename << endl;
}

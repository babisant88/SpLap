/*
 * aux.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: babis
 */

#include <iostream>
#include <fstream>

#include "../hpp/aux.hpp"

using namespace std;

void printMatrixXd(MatrixXd& X)
{
	for( int row_i=0; row_i<X.rows(); row_i++ )
	{
		cout << endl;
		for( int col_i=0; col_i<X.cols(); col_i++ )
			cout << X(row_i, col_i) << '\t';
	}
	cout << endl << endl;
}



void printVectorXd(VectorXd& X)
{
	for( int row_i=0; row_i<X.rows(); row_i++ )
		cout << X(row_i) << '\t';

	cout << endl << endl;
}



void printSpMat(SpMat& X)
{
	for(int k=0; k<X.outerSize(); ++k)
		for(SpMat::InnerIterator it(X, k); it; ++it)
			cout << "i = " << it.row() << ", j = " << it.col() << ", x = " << it.value() << endl;

	cout << endl << endl;
}



void printVectorXdToFile(VectorXd& X, string filename)
{
	ofstream myfile;
	myfile.open (filename.c_str());

	for( int row_i=0; row_i<X.rows(); row_i++ )
		myfile << X(row_i) << '\t';

	myfile.close();
}



void printMatrixXdToFile(MatrixXd& X, string filename)
{
	ofstream myfile;
	myfile.open (filename.c_str());

	for( int row_i=0; row_i<X.rows(); row_i++ )
	{
		for( int col_i=0; col_i<X.cols(); col_i++ )
			myfile << X(row_i, col_i) << '\t';

		myfile << endl;
	}

	myfile.close();
}



void printIntArrToFile(int* X, int n, string filename)
{
	ofstream myfile;
	myfile.open (filename.c_str());

	for( int i=0; i<n; i++ )
		myfile << X[i] << '\t';


	myfile.close();
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

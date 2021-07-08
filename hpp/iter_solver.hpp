/*
 * solver.hpp
 *
 *  Created on: Jul 8, 2021
 *      Author: babis
 */

#ifndef HPP_ITER_SOLVER_HPP_
#define HPP_ITER_SOLVER_HPP_

#include "types.hpp"

template<typename T>
int my_cg( T& A, VectorXd& b, VectorXd& x, double itol, int maxiter )
{
	int rho1;

	double alpha;

	VectorXd p, q;

	VectorXd r = b - A*x;

	/* Obtain the norm2 of b & r vectors */
	double r_norm = r.norm();
	double b_norm = b.norm();

	/* Set b_norm = 1 in case it's zero to avoid seg fault */
	b_norm = (b_norm == 0.0) ? 1.0 : b_norm;

	int iter = 0;

	while ( iter < maxiter && (r_norm / b_norm) > itol)
	{
		iter++;

		/* rho = r'*r */
		double rho = r.dot(r);

		if (iter == 1) {
			/* Set p = r */
			p = r;
		}
		else {
			double beta = rho / rho1;
			/* p = z + beta*p */
			p = r + beta*p;
		}

		rho1 = rho;

		q = A*p;

		/* a = rho / (p'*q) */
		alpha = rho / p.dot(q);

		/* x = x + alpha*p */
		x = x + alpha * p;

		/* r = r - alpha*q */
		r = r - alpha*q;

		r_norm = r.norm();
	}

	return iter;
}

#endif /* HPP_ITER_SOLVER_HPP_ */

// # @author: Bogdan Hlevca 2012
#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/************************************************************/
/***** Gauss-siedel method to solve system of equations *****/
/************************************************************/

using namespace boost::numeric::ublas;

struct gsResult {
	vector<double> x;
	int            kiter;
};


struct gsResult* gauss_seidel(matrix<double> A,
		 	 	   vector<double>  b,
		 	 	   vector<double>  x0,
		 	 	   int      omega,
		 	 	   int      kmax,
		 	 	   float    tol = 1.0e-9);

#endif

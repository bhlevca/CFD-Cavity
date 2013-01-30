// @author: Bogdan Hlevca 2012
#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/************************************************************/
/***** Gauss-siedel method to solve system of equations *****/
/************************************************************/

using namespace boost::numeric::ublas;

/**
 * Gauss elimination with pivoting
 * A*x = b
 */
vector<double> gauss_elimination_pivot(matrix<double> A, vector<double>  b);

/**
 * Gauss elimination w/wo pivoting
 * A*x = b
 *
 * if pivotflag = true is using pivoting , simple if pivotflag=false
 */
void  gauss_elimination(matrix<double> A, vector<double>&  b, bool pivotflag);

#endif


// # @author: Bogdan Hlevca 2012
#ifndef THOMAS_H
#define THOMAS_H

#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

/*
 * N - the size of the matrix         - in
 * b - subdiagonal vector             - in
 * a - the diagonal vector            - in
 * c - the supradiagonal vector       - in
 * x - vecor holding the solution     - out
 * q - the right hand side            - in
 *
 * return is done in x which is a reference & -in C++ reference allow in/out
 */
void ThomasAlgorithm(int N, vector<double>b, vector<double>a, vector<double>c, vector<double>&  x, vector<double>q);

#endif

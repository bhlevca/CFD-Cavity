// # @author: Bogdan Hlevca 2012
#include <stdlib.h>
#include "thomas.h"

/*
 * N - the size of the matrix         - in
 * b - subdiagonal vector             - in
 * a - the diagonal vector            - in
 * c - the supradiagonal vector       - in
 * x - vector holding the solution     - out
 * q - the right hand side            - in
 *
 * return is done in x which is a reference & -in C++ reference allow in/out
 */
void ThomasAlgorithm(int N, vector<double>b, vector<double>a, vector<double>c, vector<double>&  x, vector<double>q){
	int i;
	static vector<double> l(N);
	static vector<double> u(N);
	static vector<double> d(N);
	static vector<double> y(N);

	/* LU Decomposition */
	d(0) = a(0);
	u(0) = c(0);

	for(i=0;i<N-2;i++){
		l(i) = b(i)/d(i);
		d(i+1) = a(i+1) - l(i)*u(i);
		u(i+1) = c(i+1);
	}
	l(N-2) = b(N-2)/d(N-2);
	d(N-1) = a(N-1) - l(N-2)*u(N-2);

	/* Forward Substitution [L][y] = [q] */
	y(0) = q(0);
	for(i=1;i<N;i++)
		y(i) = q(i) - l(i-1)*y(i-1);

	/* Backward Substitution [U][x] = [y] */
	x(N-1) = y(N-1)/d(N-1);
	for(i=N-2;i>=0;i--)
		x(i) = (y(i) - u(i)*x(i+1))/d(i);

	return;
}

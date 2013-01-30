// @author: Bogdan Hlevca 2012
#include "gauss_elimination.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <boost/numeric/ublas/matrix_proxy.hpp>

static double EPSILON = 1e-10;



// Gaussian elimination with partial pivoting

vector<double> gauss_elimination_pivot(matrix<double> A, vector<double>  b) {

	int N  = b.size();

	for (int p = 0; p < N; p++) {

		// find pivot row and swap
		int max = p;
		for (int i = p + 1; i < N; i++) {
			if (fabs(A(i,p)) > fabs(A(max,p))) {
				max = i;
			}
		}
		matrix_row<matrix<double> > row_p(A, p);
		matrix_row<matrix<double> > row_max(A, max);
		row_p.swap(row_max);
		double   t = b(p); b(p) = b(max); b(max) = t;

		// check if singular or nearly singular
		if (fabs(A(p,p)) <= EPSILON) {
			throw new std::string("Matrix is singular or nearly singular");
		}

		// pivot within A and b
		for (int i = p + 1; i < N; i++) {
			double alpha = A(i,p) / A(p,p);
			b(i) -= alpha * b(p);
			for (int j = p; j < N; j++) {
				A(i,j) -= alpha * A(p,j);
			}
		}
	}

	// back substitution
	vector <double>* y = new  vector <double>(N);
	vector <double>&x =*y;
	for (int i = N - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < N; j++) {
			sum += A(i,j) * x(j);
		}
		x(i) = (b(i) - sum) / A(i,i);
	}
	return x;
}


void gauss_elimination(matrix<double> A, vector<double>&  b, bool pivotflag){
	int pivot;
	int N =  b.size();
	/* NOTE: The values contained in both the matrix A and
	the vector b are modified in this routine. Upon
	returning, A contains the upper triangular matrix
	obtained from its LU decomposition, and b contains
	the solution of the system Ax=b*/
	// Steps (1) and (2) (decomposition and solution of Ly = b)
	switch(pivotflag){
		case true: // Case in which pivoting is employed
			for(int k=0; k < N; k++){
				// find pivot row and swap
				pivot = k;
				for (int i = k + 1; i < N-1; i++) {
					if (fabs(A(i,k)) > fabs(A(pivot,k))) {
						pivot = i;
					}
				}

				matrix_row<matrix<double> > row_k(A, k);
				matrix_row<matrix<double> > row_pivot(A, pivot);
				row_k.swap(row_pivot);
				double   t = b(k); b(k) = b(pivot); b(pivot) = t;

				for(int i=k+1;i<N;i++){
					double l_ik = A(i,k)/A(k,k);
					for(int j=k;j<N;j++) {
						A(i,j) = A(i,j) - l_ik*A(k,j);
					}
					b(i) = b(i) - l_ik*b(k);
				}
			}
			break;
		case false: // Case 0/default in which no pivoting is used
		default:
			for(int k=0; k < N-1; k++){
				for(int i=k+1;i<N;i++){
					double l_ik = A(i,k)/A(k,k);
					for(int j=k;j<N;j++)
						A(i,j) = A(i,j) - l_ik*A(k,j);
					b(i) = b(i) - l_ik*b(k);
				}
			}
			break;
	}
	// Step (3) (backsolving to solve Ux=y)
	b(N-1) = b(N-1)/A(N-1,N-1);
	for(int k=N-2;k>=0;k--){
		for(int j=k+1;j<N;j++) {
			b(k) -= A(k,j)*b(j);
		}
		b(k) = b(k)/A(k,k);
	}
}

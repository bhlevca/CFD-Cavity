// @author: Bogdan Hlevca 2012

#include "gauss_seidel.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _USE_MATH_DEFINES

static double EPSILON = 1e-10;

/************************************************************/
/***** Gauss-siedel method to solve system of equations *****/
/************************************************************/

//Make the suystem diagonal dominant if possible
void pivoting_system(matrix<double>& A, vector<double>&b, int n)
{
    //pivoting
    for (int i=0; i<n; i++) {
		int pivot = i;
		for (int k = i + 1; k < n; k++) {
			if (fabs(A(k,i)) > fabs(A(pivot,i))) {
				pivot = k;
			}
		}
		matrix_row<matrix<double> > row_i(A, i);
		matrix_row<matrix<double> > row_pivot(A, pivot);
		row_i.swap(row_pivot);
		double   t = b(i); b(i) = b(pivot); b(pivot) = t;
		// singular or nearly singular
		if (fabs(A(i,i)) <= EPSILON) {
			throw new std::string("Matrix is singular or nearly singular");
		}
    }
}

struct gsResult* gauss_seidel(matrix<double> A,
		 	 	   vector<double>  b,
		 	 	   vector<double>  x0,
		 	 	   int      omega,
		 	 	   int      kmax,
		 	 	   float    tol) {


		struct  gsResult* gr = new gsResult();
		double sum = 0.0;
	 	int n = b.size();
	 	vector<double> x(x0);

	 	try {
			pivoting_system(A, b ,n);

			//iterations until maximum number of iterations or the precision is met
			for (int k=0; k< kmax; k++) {
				//copy the old values
				vector<double>xOld = x;

				//loop through all the equations
				for (int i=0; i<n; i++) {

					sum = 0;
					//loop  from  j=0 to j=i-1
					for (int j=0; j < i; j++ ) {
						sum += A(i, j) * x(j);
					}//end j loop

					//loop  from  j=i+1 to j=n-1
					for (int j = i+1; j<n; j++) {
						sum +=  A(i, j) * x(j);
					}//end j loop
					x(i) = (1.0 - omega) * x(i) + (b(i) - sum) * omega / A(i, i);
				} //end i loop

				//test the convergence error
				double dx = sqrt( inner_prod(x - xOld, x - xOld));
				if (dx < tol) {
					gr->x = x;
					gr->kiter =k;
					return gr;
				}

				// re compute relaxation coeficient omega
				//from "Parallel Scientific Computing in C++ and MPI"
				//by    George Em Karniadakis and Robert M. Kirby II
				if (k % 10 == 0) {
					omega = 2 / (1 + sin(M_PI/n));
				}
			}//end k loop
			std::cout << "Gauss-Seidel failed to converge: too many iterations" << std::endl;

	 	} catch (std::string& err) {
	 		std::cout << "Gauss-Seidel failed to converge:" << err << std::endl;
	 	}
		gr->x = x;
		gr->kiter =kmax;
		return gr;
 }


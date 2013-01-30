//# @author: Bogdan Hlevca 2012
#include "gauss_seidel.h"
#include "gauss_elimination.h"
#include "thomas.h"
#include "storage_adaptors.hpp"

#define ACC 0.00000001

int main(){

	//define the amatrix
	double initialValues [3][3] = {
			{8.0, 1.0, 6.0},
			{3.0, 5.0, 7.0},
			{4.0, 9.0, 2.0},
	};

	//define the initial values
	double initialX[3] = {1.0, 0.0, 1.0};

	//define the right hand side
	double initialB[3] = {1.0, 2.0, 3.0};


	//allocate the matrix and the vectors
	matrix<double> A(3, 3);
	vector<double> b (3);
	vector<double> x0 (3);


	//fill the values into the structures
	A = make_matrix_from_pointer(initialValues);
	b = make_vector_from_pointer(3, initialB);
	x0 = make_vector_from_pointer(3, initialX);

	//define the relaxation facto and the maximum allowed iterations
	double omega=1.0;
	int kmax=200;

	// call gauss-seidel algorithm
	struct gsResult* r=gauss_seidel(A,b,x0,omega,kmax,ACC);

	//display the GS result
	std::cout << "Gauss Seidel Iterations\n-----------------------------------------" << std::endl;
	std::cout << "x=" << r->x << std::endl;
	std::cout << "Number of iterations=" << r->kiter << std::endl;

	//re-initialize the matrix and the free term
	A = make_matrix_from_pointer(initialValues);
	b = make_vector_from_pointer(3, initialB);

	// call the Gauss elimination code
	std::cout << "\n\nGauss Elimination with partial pivoting\n------------------------------------" << std::endl;
	vector<double> x =  gauss_elimination_pivot(A,b);
	std::cout << "x=" << x << std::endl;

	std::cout << "\n\nGauss Elimination with no pivoting\n------------------------------------" << std::endl;
	gauss_elimination(A,b, false);
	std::cout << "x=" << b << std::endl;



	//Call the  TDMA algorithm
	std::cout << "\n\nThomas Algorithm\n------------------------------------" << std::endl;

	// initialize the Thomas subdiagonal vector
	double initial_B[3] = {2., 1., 1.};
	// initialize the main diagonal
	double initial_A[4] = {3., -3., 2., -1.};
	// initialize the supradiagonal vector
	double initial_C[3] = {-1., 2., 5.};
	//initialize the free terms
	double initial_Q[4] = {5, 5, 10, 1};

	//allocate space and initialize them
	vector<double> aa(4);
	aa = make_vector_from_pointer(4,initial_A);

	vector<double> bb(3);
	bb= make_vector_from_pointer(3,initial_B);

	vector<double> cc(3);
	cc= make_vector_from_pointer(3, initial_C);

	vector<double> qq(4);
	qq= make_vector_from_pointer(4,initial_Q);

	vector<double> xx(4);

	//call the algorithm
	ThomasAlgorithm(4, bb, aa, cc, xx, qq);

	//print the results
	std::cout << "x=" << xx << std::endl;
	return 0;
}


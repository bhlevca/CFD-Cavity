#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <boost/numeric/ublas/io.hpp>

double dot( int N, double* a, double* b){
	double runningSum = 0;
	for (int i=0; i< N; i++) {
		runningSum += a[i] * b[i];
	}
	std::cout << runningSum << std::endl;
	return runningSum;
}



void SOR(double omega, int N, double A[][3], double* x, double* b, double abstol){

	int maxit =1000;
	double sum1, sum2;
	double * xold = new double[N];

	//set initial guess
	for (int i=0; i<N; i++) {
		x[i] = 1.0;
		xold[i] = 1.0;
	}
	for (int k=0; k < maxit; k++) {
		for (int i=0; i<N ; i++){
			sum1=0.0; sum2=0.0;
			for(int j=i+1; j<N; j++) {
				sum2=sum2+A[i][j]*x[j];
			}
			for(int j=0;j<i; j++) {
				sum1=sum1+A[i][j]*xold[j];
			}
			x[i] = (1.0-omega)*xold[i]+omega*(b[i]-sum1-sum2)/A[i][i];
		}
		if(sqrt(dot(N,x,xold))<abstol){
			delete [] xold;
			return;
		}
		for (int i=0; i<N; i++) {
			xold[i] = x[i];
		}
	}
	std::cout << "SOR : Maximum number of iterations reaches\n" << std::endl;
	delete [] xold;
	return;
}

/*
int main()
{

	double a[3][3] = {
				{8.0, 1.0, 6.0},
				{3.0, 5.0, 7.0},
				{4.0, 9.0, 2.0}
		};
	double x[3] = {0.0, 0.0, 0.0};
	double b[3] = {1.0, 2.0, 3.0};

	SOR (0.3, 3, a, x, b, 1.0e-6);
	for (int i=0; i< 3; i++)
		std::cout << x[i] << std::endl;

}
*/

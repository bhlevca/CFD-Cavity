# @author: Bogdan Hlevca 2012
'''
    x,numIter,omega = gaussSeidel(A, b, x0, omega, kmax, tol = 1.0e-9)
    Gauss-Seidel method for solving [A]{x} = {b}.
    The matrix [A] can be sparse. 
'''

import math
import numpy

#Make the suystem diagonal dominant if possible
def pivoting_system(A, b, n):
    EPSILON = 1e-10

    #pivoting
    for i in range(0, n):
        pivot = i
        for k in range(i + 1, n):
            if math.fabs(A[k, i]) > math.fabs(A[pivot, i]):
                pivot = k
            #end if
        #end for
        if i != pivot:
            row_i = A[i, ].copy()
            row_pivot = A[pivot, ].copy()
            temp = row_i.copy()
            A[i, ] = row_pivot.copy()
            A[pivot, ] = row_i.copy()
            temp = b[i]; b[i] = b[pivot]; b[pivot] = temp

        # singular or nearly singular
        if math.fabs(A[i, i]) <= EPSILON:
            raise Exception("Matrix is singular or nearly singular")

def gaussSeidel(A, b, x0, omega, kmax, tol = 1.0e-9):
    '''
       arguments:
        A - matrix coefficients
        b - free term
        x0 - initial values
        omega - relaxation factor
        kmax - maximum number de iterations
        tol - error tolerance
        
       returns: 
           The solution vector and the number of iterations used. 
    '''
    #get the sized of the system
    n = numpy.size(b)
    x = x0

    #pivot the system if neessary
    pivoting_system(A, b , n)

    #iterations until maximum number of iterations or the precision is met
    for k in range(0, kmax):
        #copy the old values 
        xOld = x.copy()

        #loop through all the equations
        for i in range(0, n):
            sum = 0.0

            #loop  from  j=0 to j=i-1
            for j in range(0, i):
                sum = sum + A[i, j] * x[j]
                #print "i:%d, j:%d, A(i,j):%f, x[j]:%f" % (i, j, A[i, j], x[j])
            #end j loop

            #loop  from  j=i+1 to j=n-1
            for j in range(i + 1, n):
                sum = sum + A[i, j] * x[j]
                #print "i:%d, j:%d, A(i,j):%f, x[j]:%f" % (i, j, A[i, j], x[j])
            #end j loop
            x[i] = (1 - omega) * x[i] + (b[i] - sum) * omega / A[i, i]
            #print "i:%d, sum: %f, x(i): %f" % (i, sum, x[i])
        # test the convergence error 
        dx = math.sqrt(numpy.dot(x - xOld, x - xOld))
        if dx < tol:
            return x, k

    #end i loop

    print 'Gauss-Seidel failed to converge'
    return x, k




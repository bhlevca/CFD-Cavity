# @author: Bogdan Hlevca 2012 bogdan hlevca 995151213
import math
import numpy
'''
 * Algrithm from http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
 * N - the size of the matrix         - in
 * b - subdiagonal vector             - in indexed from 1..n-1
 * a - the diagonal vector            - in
 * c - the supradiagonal vector       - in indexed from 0..n-2
 * x - vector holding the solution     - out
 * q - the right hand side            - in
 
'''
def thomas(N, b, a, c, q):
    '''
    % THOMAS    Solves a tridiagonal linear system
    %
    %    x = THOMAS(A,d) solves a tridiagonal linear system using the very efficient
    %    Thomas Algorith. The vector x is the returned answer.
    %
    %       A*x = d;    /  a1  c1   0   0   0   ...   0  \   / x1 \    / q1 \
    %                   |  b1  a2  c2   0   0   ...   0  |   | x2 |    | q2 |
    %                   |   0  b2  a3  c3   0   ...   0  | x | x3 | =  | q3 |
    %                   |   :   :   :   :   :    :    :  |   | x4 |    | q4 |
    %                   |   0   0   0   0 bn-2 an-1 cn-1 |   | :  |    |  : |
    %                   \   0   0   0   0   0  bn-1  an /    \ xn /    \ qn /
    %
    %   - The matrix A must be strictly diagonally dominant for a stable solution.
    %   - This algorithm solves this system on (5n-4) multiplications/divisions and
    %      (3n-3) subtractions.
    %
    %   x = THOMAS(a,b,c,d) where a is the diagonal, b is the upper diagonal, and c is 
    %       the lower diagonal of A also solves A*x = d for x. Note that a is size n 
    %       while b and c is size n-1.
    %       If size(a)=size(d)=[L C] and size(b)=size(c)=[L-1 C], THOMAS solves the C
    %       independent systems simultaneously.
    %   
    %
    %   ATTENTION : No verification is done in order to assure that A is a tridiagonal matrix.
    %   If this function is used with a non tridiagonal matrix it will produce wrong results.
    %
    '''

    l = numpy.zeros(N)
    u = numpy.zeros(N)
    d = numpy.zeros(N)
    y = numpy.zeros(N)
    x = numpy.zeros(N)
    # LU Decomposition 
    d[0] = a[0]
    u[0] = c[0]

    for i in range(0, N - 2):
        l[i] = b[i] / d[i]
        d[i + 1] = a[i + 1] - l[i] * u[i]
        u[i + 1] = c[i + 1]
    #end for
    l[N - 2] = b[N - 2] / d[N - 2]
    d[N - 1] = a[N - 1] - l[N - 2] * u[N - 2]

    '''
    %1. LU decomposition ________________________________________________________________________________
    %
    % L = / 1                \     U =  / m1  r1              \
    %     | l1 1             |          |     m2 r2           |
    %     |    l2 1          |          |        m3 r3        |
    %     |     : : :        |          |         :  :  :     |
    %     \           ln-1 1 /          \                  mn /
    %
    %  ri = bi -> not necessary 
    '''
    #Forward Substitution [L][y] = [q] 
    y[0] = q[0]
    for i in range (1, N):
        y[i] = q[i] - l[i - 1] * y[i - 1]

    # Backward Substitution [U][x] = [y] 
    x[N - 1] = y[N - 1] / d[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = (y[i] - u[i] * x[i + 1]) / d[i]

    return x


#end

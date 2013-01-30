# @author: Bogdan Hlevca 2012

from numpy import zeros, array, prod, diagonal, dot
from gaussElimin import *
from gaussSeidel import *
from thomas import *
import timeit

# Gauss Elimination test
print "Gauss Elimination:"
print "------------------"
print "A =              b = "
print "8.0, 1.0, 6.0      1 "
print "3.0, 5.0, 7.0      2 "
print "4.0, 9.0, 2.0      3 "

a = zeros((3, 3))
a[0] = array([8.0, 1.0, 6.0])
a[1] = array([3.0, 5.0, 7.0])
a[2] = array([4.0, 9.0, 2.0])

b = array([1.0, 2.0, 3.0])

aOrig = a.copy() # Save original matrix
bOrig = b.copy() # and the constant vector

x = gaussElimin(a, b)

det = prod(diagonal(a))
print 'Solution: x =\n', x
print '\ndet =', det
print '\nCheck result: [a]{x} - b =\n', dot(aOrig, x) - bOrig

print "\n\n-----------------------"
print "Gauss Seidel:"
print "------------------"
print "A =              b = "
print "8.0, 1.0, 6.0      1 "
print "3.0, 5.0, 7.0      2 "
print "4.0, 9.0, 2.0      3 "
#Gauss Seidel test
t = timeit.Timer("print 'initialize + solve:'", "print 'Python Gauss-Seidel'")

a = numpy.zeros((3, 3))
a[0] = numpy.array([8.0, 1.0, 6.0])
a[1] = numpy.array([3.0, 5.0, 7.0])
a[2] = numpy.array([4.0, 9.0, 2.0])

b = numpy.array([1.0, 2.0, 3.0])

#a[0] = numpy.array([12.0, 3.0, -5.0])
#a[1] = numpy.array([1.0, 5.0, 3.0])
#a[2] = numpy.array([3.0, 7.0, 13.0])

#b = numpy.array([1.0, 28.0, 76.0])

x0 = numpy.array([1.0, 0.0, 1.0])

y, kiter = gaussSeidel(a, b, x0, 1.0, 20)

print "Timer:", t.timeit(1)
print "\nSolution:", y
print "-----------------------"
print "Number of iterations: %d" % kiter


print "\n\n-----------------------"
print "Gauss Elim & Gauss Seidel:"
print '''
    a[0] = numpy.array([20, -5, 0, 0, 0])
    a[1] = numpy.array([-5, 15, -5, 0, 0])
    a[2] = numpy.array([0, -5, 15, -5, 0])
    a[3] = numpy.array([0, 0, -5, 15, -5])
    a[4] = numpy.array([0, 0, 0, -5, 10])
    b = numpy.array([1100, 100.0, 100.0, 100, 100])
    '''
a = zeros((5, 5))
a[0] = numpy.array([20, -5, 0, 0, 0])
a[1] = numpy.array([-5, 15, -5, 0, 0])
a[2] = numpy.array([0, -5, 15, -5, 0])
a[3] = numpy.array([0, 0, -5, 15, -5])
a[4] = numpy.array([0, 0, 0, -5, 10])
b = numpy.array([1100, 100.0, 100.0, 100, 100])



aOrig = a.copy() # Save original matrix
bOrig = b.copy() # and the constant vector
x = gaussElimin(a, b)

det = prod(diagonal(a))
print 'Gauss Elim: x =\n', x
print '\ndet =', det
print '\nCheck result: [a]{x} - b =\n', dot(aOrig, x) - bOrig


x0 = numpy.array([1.0, 0.0, 1.0, 0.0, 1.0])
y, kiter = gaussSeidel(aOrig, bOrig, x0, 1.0, 500)


print "\n-----------------------"
print "Gauss Seidel Solution:", y
print "-----------------------"
print "Number of iterations: %d" % kiter


#Tridiagonal
print ''' \n\n Tridiagonal:
----------------------
     3 -1  0  0       5
A =  2 -3  2  0  b =  5
     0  1  2  5      10
     0  0  1 -1       1

a2,a3, a4 = 2,1,1
b1,b2,b3,b4 = 3,-3,2,-1
c1,c2,c3= -1,2,5

solution:2,1,2,1
'''
t = timeit.Timer("print 'initialize + solve:'", "print 'Python TDMA QUICK'")
bb = (2., 1., 1.)
aa = (3., -3., 2., -1.)
cc = (-1., 2., 5.)
dd = (5, 5, 10, 1)

bb = (-0.7, -0.675, -0.675, -0.817)
aa = (2.175, 1.075, 1.075, 1.075, 1.925)
cc = (-0.592, -0.425, -0.425, -0.425)
dd = (1.583, -0.05, 0, 0, 0)

a = numpy.array(aa)
b = numpy.array(bb)
c = numpy.array(cc)
d = numpy.array(dd)



x = thomas(a.size, b, a, c, d)
print "Timer:", t.timeit(1)
print "Solution:",
print x

print "\n\n\n"
t = timeit.Timer("print 'initialize + solve:'", "print 'Python Gauss-Seidel'")
a = zeros((5, 5))
TA = 100
TB = 500
a[0] = numpy.array([300, -100, 0, 0, 0])
a[1] = numpy.array([-100, 200, -100, 0, 0])
a[2] = numpy.array([0, -100, 200, -100, 0])
a[3] = numpy.array([0, 0, -100, 200, -100])
a[4] = numpy.array([0, 0, 0, -100, 300])
b = numpy.array([200 * TA, 0, 0, 0, 200 * TB])
x = numpy.array([0, 0, 0, 0, 0])
y, kiter = gaussSeidel(a, b, x0, 1.0, 200)

print "Timer:", t.timeit(1)
print "\nSolution:", y
print "-----------------------"
print "Number of iterations: %d" % kiter

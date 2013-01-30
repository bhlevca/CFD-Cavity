# @author: Bogdan Hlevca 2012
import numpy
import linalg
import timeit

print "Linalg"
print "---------"

# Create a mesh class that holds a vector of nodes
cls = linalg.CppBlas()

t = timeit.Timer("print 'initialize + solve:'", "print 'C++ Gauss-Seidel'")
a = numpy.zeros((3, 3))
a[0] = numpy.array([8.0, 1.0, 6.0])
a[1] = numpy.array([3.0, 5.0, 7.0])
a[2] = numpy.array([4.0, 9.0, 2.0])

cls.setA(a, 3, 3)
b = numpy.array([1.0, 2.0, 3.0])
cls.setB(b, b.size)

x0 = numpy.array([1.0, 0.0, 1.0])
cls.setX0(x0, x0.size)

cls.setParams(1.0, 200, 0.00000001)

x = numpy.zeros(b.size)
d = cls.solve(x, b.size)

print "Timer:", t.timeit(1)
print "Solution:", d["solution"]
print "Iterations:", d["kiter"]


t = timeit.Timer("print 'initialize + solve:'", "print 'Python Gauss-Seidel'")
a = numpy.zeros((5, 5))
TA = 100.0
TB = 500.0
a[0] = numpy.array([300, -100, 0, 0, 0])
a[1] = numpy.array([-100, 200, -100, 0, 0])
a[2] = numpy.array([0, -100, 200, -100, 0])
a[3] = numpy.array([0, 0, -100, 200, -100])
a[4] = numpy.array([0, 0, 0, -100, 300])
b = numpy.array([200 * TA, 0., 0., 0., 200 * TB])
x0 = numpy.array([0., 0., 0., 0., 0.])

cls.setA(a, 5, 5)
cls.setB(b, b.size)
cls.setX0(x0, x0.size)
cls.setParams(1.0, 200, 0.000000001)

x = numpy.zeros(b.size)
d = cls.solve(x, b.size)

print "Timer:", t.timeit(1)
print "Solution:", d["solution"]
print "Iterations:", d["kiter"]


print
print "TDMA from C++"

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
t = timeit.Timer("print 'initialize + solve:'", "print 'C++ TDMA'")
bb = (2., 1., 1.)
aa = (3., -3., 2., -1.)
cc = (-1., 2., 5.)
dd = (5.0, 5.0, 10.0, 1)

a = numpy.array(aa)
b = numpy.array(bb)
c = numpy.array(cc)
d = numpy.array(dd)
n = a.size

x = numpy.zeros(n)

cls.setTDMA(b, a, c, d, n)
d = cls.solveTDMA(x, n)
print "Timer:", t.timeit(1)
print "Solution:", d["solution"]



print
print "TDMA from C++"

print ''' \n\n Tridiagonal from QUICK:
----------------------
'''
t = timeit.Timer("print 'initialize + solve:'", "print 'C++ TDMA'")
bb = (-0.7, -0.675, -0.675, -0.817)
aa = (2.175, 1.075, 1.075, 1.075, 1.925)
cc = (-0.592, -0.425, -0.425, -0.425)
dd = (1.583, -0.05, 0, 0, 0)

a = numpy.array(aa)
b = numpy.array(bb)
c = numpy.array(cc)
d = numpy.array(dd)
n = a.size

x = numpy.zeros(n)

cls.setTDMA(b, a, c, d, n)
d = cls.solveTDMA(x, n)
print "Timer:", t.timeit(1)
print "Solution:", d["solution"]


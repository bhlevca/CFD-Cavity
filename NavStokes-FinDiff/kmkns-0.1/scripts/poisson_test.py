#!/usr/bin/env python

from numpy import *
from kmkns.poisson import *

print "Initializing arrays..."

#n = (200, 200) # with 1GB of memory, about 500k nodes is managable.
n = (20, 30, 40)

# set up some arbitrary 'f'
f = array(range(n[0] * n[1] * n[2]), float)
f.shape = n

delta = (0.001, 0.002, 0.003)
mask = ones(n, int)
# make solution unique
remove_singularity(mask)
# make equation solvable
f[1:-1,1:-1,1:-1] -= f[1:-1,1:-1,1:-1].mean()

# set up some arbitrary 'rho'
rho = f.copy() % 13 + 1

phi = zeros(n,float)

print "Solving Poisson equation..."

poisson(phi, f, delta, mask, rho)

print "Calculating residual..."

rx = (2.0 / delta[0]**2) * \
    ((phi[2:,1:-1,1:-1] - phi[1:-1,1:-1,1:-1]) / (rho[2:,1:-1,1:-1] + rho[1:-1,1:-1,1:-1]) - \
    (phi[1:-1,1:-1,1:-1] - phi[:-2,1:-1,1:-1]) / (rho[1:-1,1:-1,1:-1] + rho[:-2,1:-1,1:-1])) 

ry = (2.0 / delta[1]**2) * \
    ((phi[1:-1,2:,1:-1] - phi[1:-1,1:-1,1:-1]) / (rho[1:-1,2:,1:-1] + rho[1:-1,1:-1,1:-1]) - \
    (phi[1:-1,1:-1,1:-1] - phi[1:-1,:-2,1:-1]) / (rho[1:-1,1:-1,1:-1] + rho[1:-1,:-2,1:-1])) 

rz = (2.0 / delta[2]**2) * \
    ((phi[1:-1,1:-1,2:] - phi[1:-1,1:-1,1:-1]) / (rho[1:-1,1:-1,2:] + rho[1:-1,1:-1,1:-1]) - \
    (phi[1:-1,1:-1,1:-1] - phi[1:-1,1:-1,:-2]) / (rho[1:-1,1:-1,1:-1] + rho[1:-1,1:-1,:-2])) 

print "Done."
print "f.shape = %s" % str(f.shape)
print "rx.shape = %s" % str(rx.shape)
print "ry.shape = %s" % str(ry.shape)
print "rz.shape = %s" % str(rz.shape)
print "L2-norm of f:", (f[1:-1,1:-1,1:-1]**2).sum()**0.5
print "L2-norm of rx:", (rx**2).sum()**0.5
print "L2-norm of ry:", (ry**2).sum()**0.5
print "L2-norm of rz:", (rz**2).sum()**0.5
print "L2-norm of rx + ry + rz:", ((rx + ry + rz)**2).sum()**0.5
print "L2-norm of residual:", ((f[1:-1,1:-1,1:-1] + rx + ry + rz)**2).sum()**0.5


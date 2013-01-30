#!/usr/bin/env python

from kmkns.NavierStokes import *
import sys, getopt
import time

def run(size):
    t0 = time.clock()
    grid = Grid([0.0, 0.0], [2.0, 4.0], [size, 2 * size])
    ns = NavierStokes2D(grid, mu=0.002, mu_gas=3.2e-5,
        rho=1.0, rho_gas=0.0013, surface_tension_coeff=1.46, gravity=[0.0, -5.0])

    ns.set_bc('west+east+south+north', u=0, v=0)
    ns.add_bubble([1.0, 1.0], 0.5)

    while ns.t < 0.1:
        ns.step(ns.find_suitable_dt(), beta=0.1337)
    t1 = time.clock()
    print "Total time: %gs" % (t1 - t0)


opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
size = 50
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = float(val)

run(size)

import cProfile
cProfile.run("run(%d)" % size, sort=2)
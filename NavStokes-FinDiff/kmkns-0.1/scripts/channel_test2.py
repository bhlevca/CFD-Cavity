#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'r:', ['Re='])
filename = "channel_test2"
Re = 1.0
for opt, val in opts:
    if opt in ('-r', '--Re'):
        Re = float(val)

sizes = [6, 9, 12, 18, 25, 37, 50]

kmkns.Keyboard.set_poll()

for size in sizes:
    grid = Grid([0.0, 0.0], [3.0, 1.0], [3 * size, 2 * size + 1])
    ns = NavierStokes2D(grid, nu=1.0/Re)

    ns.set_bc('north+south', u=0, v=0)
    ns.set_bc('west', u=1, dvdn=0)
    ns.set_bc('east', dvdn=0, dudn=0)

    # try to load file
    ns.load("%s_%i.dat" % (filename, size))

    def show_info():
        print "Iteration #%i" % ns.it
        print "Time: %f" % ns.t
        mid = ns.u.shape[1]/2
        print "Centre velocitiy (u):", ns.u[-1,mid]
        print "Error: %f" % (1.5 - ns.u[-1,mid])
        ns.plot_velocity(1)

    while ns.t < 0.25:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break

        dt = min(ns.find_suitable_dt(), 0.02)
        ns.step(dt)

        if ch == 'p' or ch == 'P':
            show_info()

        if (ch == 'l' or ch == 'L'):
            ns.load("%s_%i.dat" % (filename, size))
        if (ch == 's' or ch == 'S'):
            ns.save("%s_%i.dat" % (filename, size))

    if ns.t < 0.25:
        break

    print "Done with %ix%i." % (grid.shape[0] - 1, grid.shape[1] - 1)
    show_info()
    
    ns.save("%s_%i.dat" % (filename, size))

kmkns.Keyboard.set_normal()

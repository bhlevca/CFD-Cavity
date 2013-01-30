#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

#opts, args = getopt.getopt(sys.argv[1:], 'r:', ['Re='])
filename = "advection_test2"
#Re = 1.0
#for opt, val in opts:
#    if opt in ('-r', '--Re'):
#        Re = float(val)

sizes = [50, 100, 150, 200]

kmkns.Keyboard.set_poll()

for size in sizes:
    grid = Grid([-1.0, -1.0], [1.0, 1.0], [size, size])
    ns = NavierStokes2D(grid, mu_liquid=1.0, mu_gas=1.0, rho_liquid=1.0, rho_gas=1.0)

    u_bc = lambda x, y, t: -y
    v_bc = lambda x, y, t: x
    u_ic = lambda x, y: -y
    v_ic = lambda x, y: x
    ns.set_ic(u=u_ic, v=v_ic)
    ns.set_bc('north+south+west+east', u=u_bc, v=v_bc)
    
    ns.add_bubble([0, 0.5], 0.3)
    ns.add_liquid([-0.05, -1.0], [0.05, 0.7])

    def show_info():
        print "Iteration #%i" % ns.it
        print "Time: %f" % ns.t
        ns.plot_level([0.5], 1)

    show_info()
    if viz.backend == 'gnuplot':
        g = viz.get_backend()
        g("set nokey")
        g("set size square")
        g("set xrange [-1:1]")
        g("set yrange [-1:1]")
        g("set rmargin 12")
        g("set lmargin 12")
        g("set bmargin 12")
        g("set tmargin 12")
        g("replot")
        g.hardcopy("advection%i_init.eps" % size)
    else:
        viz.hardcopy("advection%i_init.eps" % size)


    # try to load file
    ns.load("%s_%i.dat" % (filename, size))

    while ns.t < 2.0 * pi:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break

        ns.step(ns.find_suitable_dt())

        if ch == 'p' or ch == 'P':
            show_info()

        if (ch == 'l' or ch == 'L'):
            ns.load("%s_%i.dat" % (filename, size))
        if (ch == 's' or ch == 'S'):
            ns.save("%s_%i.dat" % (filename, size))

    if ns.t < 2.0 * pi:
        break

    print "Done with %ix%i." % (grid.shape[0] - 1, grid.shape[1] - 1)
    show_info()
    if viz.backend == 'gnuplot':
        g = viz.get_backend()
        g("set nokey")
        g("set size square")
        g("set xrange [-1:1]")
        g("set yrange [-1:1]")
        g("set rmargin 12")
        g("set lmargin 12")
        g("set bmargin 12")
        g("set tmargin 12")
        g("replot")
        g.hardcopy("advection%i_final.eps" % size)
    else:
        viz.hardcopy("advection%i_final.eps" % size)

    ns.save("%s_%i.dat" % (filename, size))

kmkns.Keyboard.set_normal()

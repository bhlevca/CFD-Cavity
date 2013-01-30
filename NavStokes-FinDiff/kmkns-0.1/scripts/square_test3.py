#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

def grid_range(t):
    if len(t) > 1:
        for i in xrange(t[0]):
            for j in grid_range(t[1:]):
                yield (i,) + j
    else:            
        for i in xrange(t[0]):
            yield (i,)

def save_plot(filename):
    if viz.backend == "gnuplot":
        g = viz.get_backend()
        g("unset key")
        g("set size square")
        g("set xrange [0:1]")
        g("set yrange [0:1]")
        g("set rmargin 12")
        g("set lmargin 12")
        g("set bmargin 12")
        g("set tmargin 12")
        g("replot")
        # viz.hardcopy() does not work with custom settings!!!
        g.hardcopy(filename)
    else:
        viz.hardcopy(filename)

opts, args = getopt.getopt(sys.argv[1:], 'p:', ['power='])
filename = "square_test3"
power = 7
for opt, val in opts:
    if opt in ('-p', '--power'):
        power = int(val)

Res = (1.0, 10.0, 100.0, 1000.0, 5000.0)
times = (0.5, 5.0, 50.0, 390.0, 1600.0)

stream_contours = [-0.1, -0.08, -0.06, -0.04, -0.02, -0.01,
    -3e-3, -1e-3, -3e-4, -1e-4, -3e-5, -1e-5, -3e-6, -1e-6,
    -1e-7, -1e-8, -1e-9, -1e-10, 0.0, 1e-10, 1e-9, 1e-8, 1e-7,
    1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03,
    0.05, 0.07, 0.09, 0.1, 0.11, 0.115, 0.1175]

vorticity_contours = [-40.0, -35.0, -30.0, -25.0, -20.0, -15.0,
    -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0, -0.5, -0.2,
    0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    15.0, 20.0, 25.0, 30.0, 35.0, 40.0]

grid = Grid([0.0, 0.0], [1.0, 1.0], [2**power, 2**power])

kmkns.Keyboard.set_poll()

for Re, time in zip(Res, times):
    ns = NavierStokes2D(grid, nu=1.0/Re)
    ns.set_bc('west+east+south', u=0, v=0)
    ns.set_bc('north', u=-1, v=0)

    # try to load file
    ns.load("%s_%g.dat" % (filename, Re))
    
    def show_info(save):
        ns.plot_stream(stream_contours, 1)
        viz.title("Stream contours, Re=%g" % Re)
        if save:
            save_plot("square%gs_inf.eps" % Re)
        
        ns.plot_vorticity(vorticity_contours, 2)
        viz.title("Vorticity contours, Re=%g" % Re)
        if save:
            save_plot("square%gv_inf.eps" % Re)

        print "Iteration #%i" % ns.it
        print "Time: %f" % ns.t
        stream = ns.stream_function()
        vort = ns.vorticity()
        smax = stream.max()
        for (i, j) in grid_range(stream.shape):
            if stream[i,j] == smax:
                print "Max stream:", smax, "at", grid[i,j]
                print "Vorticity:", vort[i,j]
                break
        smin = stream.min()
        for (i, j) in grid_range(stream.shape):
            if stream[i,j] == smin:
                print "Min stream:", smin, "at", grid[i,j]
                print "Vorticity:", vort[i,j]
                break

    while ns.t < time:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break
            
        ns.step(ns.find_suitable_dt())

        if ch == 'p' or ch == 'P':
            show_info(False)

        if (ch == 'l' or ch == 'L'):
            ns.load("%s_%g.dat" % (filename, Re))
        if (ch == 's' or ch == 'S'):
            ns.save("%s_%g.dat" % (filename, Re))

    if ns.t < time:
        break

    print "Done with Re=%g" % Re

    show_info(True)

    ns.save("%s_%g.dat" % (filename, Re))

dummy_a = 1
dummy_b = 1
for dummy_loop in xrange(200000):
    dummy_c = dummy_a + dummy_b
    dummy_a = dummy_b
    dummy_b = dummy_c

kmkns.Keyboard.set_normal()

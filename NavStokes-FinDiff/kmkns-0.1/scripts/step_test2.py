#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
import kmkns.Keyboard

def grid_range(t):
    if len(t) > 1:
        for i in xrange(t[0]):
            for j in grid_range(t[1:]):
                yield (i,) + j
    else:            
        for i in xrange(t[0]):
            yield (i,)

opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
filename = "step_test2"
size = 100
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)

Res = (100.0, 200.0, 300.0, 400.0, 500.0, 600.0)
times = (50.0, 100.0, 150.0, 200.0, 250.0, 300.0)

stream_contours = [0.003 * i for i in xrange(-50, 0)] + \
    [-0.001, -0.0003, -0.0001, 0.0, 0.0001, 0.0003, 0.001, 0.003] + \
    [0.02 * i for i in xrange(1, 25)] + \
    [0.497, 0.499, 0.4997, 0.4999, 0.5, 0.5001] + \
    [0.5 + 0.0003 * i for i in xrange(1, 50)]

grid = Grid([0.0, 0.0], [15.0, 1.0], [3 * size, size])
    
kmkns.Keyboard.set_poll()

for Re, time in zip(Res, times):
    ns = NavierStokes2D(grid, nu=1.0/Re)
    ns.set_bc('north+south', u=0, v=0)
    ns.set_bc('east', dudn=0, dvdn=0)
    ns.set_bc('west', u=lambda x, y, t: max(0.0, 24.0 * (1.0 - y) * (y - 0.5)), v=0)

    # try to load file
    ns.load("%s_%g.dat" % (filename, Re))

    def show_info(save):
        ns.plot_stream(stream_contours, 1)
        viz.title("Stream contours, Re=%g" % Re)
        if save:
            if viz.backend == "gnuplot":
                g = viz.get_backend()
                g("unset key")
                g("set size ratio 0.2")
                g("set xrange [0:15]")
                g("set yrange [0:1]")
                g("set rmargin 10")
                g("set lmargin 10")
                g("set bmargin 10")
                g("set tmargin 10")
                g("replot")
                # viz.hardcopy() does not work with custom settings!!!
                g.hardcopy("step%g.eps" % Re)
            else:
                viz.hardcopy("step%g.eps" % Re)

        print "Iteration #", ns.it
        print "Time:", ns.t
        stream = ns.stream_function()
        smax = stream.max()
        smin = stream.min()
        for (i, j) in grid_range(stream.shape):
            if smin == stream[i,j]:
                print "Min stream:", smin, "at", grid[i,j]
                break
        for (i, j) in grid_range(stream.shape):
            if smax == stream[i,j]:
                print "Max stream:", smax, "at", grid[i,j]
                break
        x = ns.p.shape[0]
        for i in xrange(1, ns.u.shape[0]):
            if ns.u[i,1] < 0.0:
                print "Lower detach./h =", 2.0 * (grid[i,0][0] - 0.5 * grid.delta[0])
                #print "Lower detach./h =", 1.5*(grid[i,0][0] - 0.5*grid.delta[0] - 7.5)
                x = i
                break
        for i in xrange(x, ns.u.shape[0]):
            if ns.u[i,1] > 0.0:
                print "Lower attach./h =", 2.0 * (grid[i,0][0] - 0.5 * grid.delta[0])
                #print "Lower attach./h =", 1.5*(grid[i,0][0] - 0.5*grid.delta[0] - 7.5)
                break
        x = ns.u.shape[0]
        for i in xrange(1, ns.u.shape[0]):
            if ns.u[i,-2] < 0.0:
                print "Upper detach./h =", 2.0 * (grid[i,0][0] - 0.5 * grid.delta[0])
                #print "Upper detach./h =", 1.5*(grid[i,0][0] - 0.5*grid.delta[0] - 7.5)
                x = i
                break
        for i in xrange(x, ns.u.shape[0]):
            if ns.u[i,-2] > 0.0:
                print "Upper attach./h =", 2.0 * (grid[i,0][0] - 0.5 * grid.delta[0])
                #print "Upper attach./h =", 1.5*(grid[i,0][0] - 0.5*grid.delta[0] - 7.5)
                break
        #ns.plot_velocity_profile(0, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0], 2)
        #viz.title("Velocity field, Re=%g" % Re)
        

    while ns.t < time:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break
            
        dt = min(ns.find_suitable_dt(), 0.5)
        ns.step(dt)

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

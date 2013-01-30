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

opts, args = getopt.getopt(sys.argv[1:], 's:l:r:w:', \
    ['save=','load=','Re=','size='])
loadfile = ""
savefile = ""
Re = 100.0
size = 50
for opt, val in opts:
    if opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val
    elif opt in ('-r', '--Re'):
        Re = float(val)
    elif opt in ('-w', '--size'):
        size = int(val)

grid = Grid([0.0, 0.0], [15.0, 1.0], [3 * size, size])
ns = NavierStokes2D(grid, nu=1.0/Re)

ns.set_bc('north+south', u=0, v=0)
ns.set_bc('east', dudn=0, dvdn=0)
ns.set_bc('west', u=lambda x, y, t: max(0.0, 24.0 * (1.0 - y) * (y - 0.5)), v=0)
ns.add_tracer([-1.0,0.5], [0.0,0.55])
ns.add_tracer([-1.0,0.6], [0.0,0.65])
ns.add_tracer([-1.0,0.7], [0.0,0.75])
ns.add_tracer([-1.0,0.8], [0.0,0.85])
ns.add_tracer([-1.0,0.9], [0.0,0.95])

kmkns.Keyboard.set_poll()
while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break

    ns.step(ns.find_suitable_dt(), beta=0.0, gamma=0.9)
    if ch == 'p' or ch == 'P':
        #ns.plot_pressure(20, 1)
        ns.plot_stream(20, 1)
        ns.plot_velocity(2)
        ns.plot_tracer(3)
        
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
    
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)
    
kmkns.Keyboard.set_normal()

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

opts, args = getopt.getopt(sys.argv[1:], 's:l:r:p:', ['save=','load=','Re=', 'power='])
loadfile = ""
savefile = ""
power = 6
Re = 1000.0
for opt, val in opts:
    if opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val
    elif opt in ('-r', '--Re'):
        Re = float(val)
    elif opt in ('-p', '--power'):
        power = int(val)

times = (1.4, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, 8.4, 9.4, 12.6, 15.0)

grid = Grid([0.0, 0.0], [1.0, 1.0], [2**power, 2**power])
ns = NavierStokes2D(grid, nu=1.0/Re)
ns.set_bc('west+east+south', u=0, v=0)
ns.set_bc('north', u=-1, v=0)

stream_contours = [-0.1, -0.08, -0.06, -0.04, -0.02, -0.01, \
    -3e-3, -1e-3, -3e-4, -1e-4, -3e-5, -1e-5, -3e-6, -1e-6, \
    -1e-7, -1e-8, -1e-9, -1e-10, 0.0, 1e-10, 1e-9, 1e-8, 1e-7, \
    1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, \
    0.05, 0.07, 0.09, 0.1, 0.11, 0.115, 0.1175]

vorticity_contours = [-40.0, -35.0, -30.0, -25.0, -20.0, -15.0, \
    -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0, -0.5, -0.2, \
    0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, \
    15.0, 20.0, 25.0, 30.0, 35.0, 40.0]


kmkns.Keyboard.set_poll()
time_index = 0
while ns.t < 17.0:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
    dt = ns.find_suitable_dt()
    if(time_index < len(times) and (ns.t + dt >= times[time_index])):
        dt = times[time_index] - ns.t
        ns.step(dt)

        ns.plot_stream(stream_contours, time_index + 1)
        viz.title("t=%g" % ns.t)
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
            g.hardcopy("square%i_%i_%i.eps" % (int(round(Re)), int(ns.t), int(ns.t*10-int(ns.t)*10)))
        
        #ns.save("square%i_%i_%i.dat" % (int(round(Re)), int(ns.t), int(ns.t*10-int(ns.t)*10)))
        time_index += 1
    else:
        ns.step(dt)
        
    if ch == 'p' or ch == 'P':
        #print ns.velocity[0,-2,:].max()
        #ns.plotVelocity(gpu)
        #ns.plotPressure(gpp)
        #ns.plotStream(gpp)
        ns.plot_stream(stream_contours, time_index + 1)
        viz.title("t=%g" % ns.t)
        
        print "Iteration #", ns.it
        print "Time:", ns.t
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
    
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)
        
kmkns.Keyboard.set_normal()

#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
filename = "drop_test2"
size = 100
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)

grid = Grid([0.0, 0.0], [6.0, 6.0], [size, size])
ns = NavierStokes2D(grid, mu=0.1, mu_gas=0.0016,
    rho=1.0, rho_gas=0.0013, surface_tension_coeff=730.0, 
    gravity=[0.0, -100.0], curv_method="mdac",
    property_smoothing=False)

ns.set_bc('west+east', u=0, v=0)
ns.set_bc('south+north', u=0, v=0)
ns.add_gas([0.0, 0.0], [6.0, 6.0])
ns.add_liquid([0.0, 0.0], [6.0, 2.0])
ns.add_drop([3.0, 4.5], 0.3)

ns.load("%s.dat" % filename)

kmkns.Keyboard.set_poll()

def plot_plots(save):
    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "# of cells > 0.5 = %d" % (ns.c > 0.5).sum()
    print "Total mass: %f" % ns.c[1:-1,1:-1].sum()

    ns.plot_velocity(fig=3)
    ns.plot_pressure(20, fig=2)
    ns.plot_level([0.5], fig=1)
    if save:
        filename = "drop_%s.eps" % str(ns.t).replace('.', '_')
        viz.title("t=%g" % ns.t)
        if viz.backend == 'gnuplot':
            g = viz.get_backend()
            g("set size square")
            g("set xrange [0:6]")
            g("set yrange [0:6]")
            g("set rmargin 23")
            g("set lmargin 23")
            g("set bmargin 23")
            g("set tmargin 23")
            g("unset key")
            #g("replot")
            g.hardcopy(filename)
        else:
            viz.hardcopy(filename)

plot_time_index = 0
plot_times = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18,
    0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.202, 0.204, 0.206,
    0.22, 0.24, 0.26, 0.28, 0.3]

while plot_time_index < len(plot_times) and plot_times[plot_time_index] <= ns.t:
    plot_time_index += 1

# plot initial state
if ns.t == 0.0:
    plot_plots(True)

while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    dt = ns.find_suitable_dt()
    if plot_time_index < len(plot_times) and \
        (ns.t + dt >= plot_times[plot_time_index]):
        ns.step(plot_times[plot_time_index] - ns.t)
        plot_plots(True)
        plot_time_index += 1
    else:
        ns.step(dt)
    
    if ch == 'p' or ch == 'P':
        plot_plots(False)
    if ch == 'l' or ch == 'L':
        ns.load("%s.dat" % filename)
    if ch == 's' or ch == 'S':
        ns.save("%s.dat" % filename)


kmkns.Keyboard.set_normal()

ns.save("%s.dat" % filename)
        

#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
filename = "rayleigh_taylor_test2"
size = 128
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)

grid = Grid([0.0, 0.0], [1.0, 4.0], [size, 4 * size])
ns = NavierStokes2D(grid, mu=0.00313, mu_gas=0.00313, \
    rho=1.225, rho_gas=0.1694, surface_tension_coeff=0.0, gravity=[0.0, -9.81])

ns.set_bc('west+east', u=0, dvdn=0)
ns.set_bc('south+north', u=0, v=0)
# f(x) = 2 - 0.05 * cos(2 * pi * x)
# F(x) = 2 * x - 0.05 * sin(2 * pi * x) / (2 * pi)
ns.set_interface(2, lambda x: 2 * x - 0.05 * sin(2 * pi * x) / (2 * pi))

ns.load("%s.dat" % filename)

level_contours = [0.5]

kmkns.Keyboard.set_poll()

def plot_plots(save):
    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "# of cells > 0.5 = %d" % (ns.c > 0.5).sum()
    print "Total mass: %f" % ns.c[1:-1,1:-1].sum()

    #ns.plot_velocity(fig=2)
    #ns.plot_pressure(20, fig=3)
    ns.plot_level(level_contours, fig=1)
    if save:
        filename = "rayleigh_taylor_%s.eps" % str(ns.t).replace('.', '_')
        viz.title("t=%g" % ns.t)
        if viz.backend == 'gnuplot':
            g = viz.get_backend()
            g("set size ratio 4")
            g("set xrange [0:1]")
            g("set yrange [0:4]")
            g("set rmargin 12")
            g("set lmargin 12")
            g("set bmargin 12")
            g("set tmargin 12")
            g("set xtic 0.5")
            g("unset key")
            g("replot")
            g.hardcopy(filename)
        else:
            viz.hardcopy(filename)

plot_time_index = 0
plot_times = [0.7, 0.8, 0.9]

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
        

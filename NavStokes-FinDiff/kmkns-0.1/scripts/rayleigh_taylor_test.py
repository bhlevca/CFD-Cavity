#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:s:l:', ['size=', 'save=', 'load='])
savefile = ""
loadfile = ""
size = 64
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)
    elif opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val

grid = Grid([0.0, 0.0], [1.0, 4.0], [size, 4 * size])
ns = NavierStokes2D(grid, mu=0.00313, mu_gas=0.00313,
    rho=1.225, rho_gas=0.1694, surface_tension_coeff=0.0, gravity=[0.0, -9.81])

ns.set_bc('west+east', u=0, dvdn=0)
ns.set_bc('south+north', u=0, v=0)
# f(x) = 2 - 0.05 * cos(2 * pi * x)
# F(x) = 2 * x - 0.05 * sin(2 * pi * x) / (2 * pi)
ns.set_interface(2, lambda x: 2 * x - 0.05 * sin(2 * pi * x) / (2 * pi))

level_contours = [0.5]

kmkns.Keyboard.set_poll()

def plot_plots(save):
    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "# of cells > 0.5 = %d" % (ns.c > 0.5).sum()
    print "Total mass: %f" % ns.c[1:-1,1:-1].sum()

    ns.plot_level(level_contours, fig=1)
    #ns.plot_velocity(fig=2)
    #ns.plot_pressure(20, fig=3)

while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    ns.step(ns.find_suitable_dt())
    
    if ch == 'p' or ch == 'P':
        plot_plots(False)
    if ch == 'l' or ch == 'L':
        ns.load("%s.dat" % filename)
    if ch == 's' or ch == 'S':
        ns.save("%s.dat" % filename)


kmkns.Keyboard.set_normal()
        

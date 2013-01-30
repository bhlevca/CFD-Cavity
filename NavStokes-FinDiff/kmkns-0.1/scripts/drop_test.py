#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'l:s:w:', ['load=', 'save=', 'size='])
loadfile = ""
savefile = ""
size = 100
for opt, val in opts:
    if opt in ('-l', '--load'):
        load = val
    elif opt in ('-s', '--save'):
        save = val
    elif opt in ('-w', '--size'):
        size = int(val)

grid = Grid([0.0, 0.0], [6.0, 6.0], [size, size])
ns = NavierStokes2D(grid, mu=0.1, mu_gas=0.0016, \
    rho=1.0, rho_gas=0.0013, surface_tension_coeff=730.0, gravity=[0.0, -100.0])

ns.set_bc('west+east', u=0, v=0)
ns.set_bc('south+north', u=0, v=0)
ns.add_gas([0.0, 0.0], [6.0, 6.0])
ns.add_liquid([0.0, 0.0], [6.0, 2.0])
ns.add_drop([3.0, 4.5], 0.3)

level_contours = [0.3, 0.5, 0.7]

kmkns.Keyboard.set_poll()

def plot_plots():
    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "# of cells > 0.5 = %d" % (ns.c > 0.5).sum()
    print "Total mass: %f" % ns.c[1:-1,1:-1].sum()

    ns.plot_level(level_contours, 1)
    ns.plot_velocity(2)
    ns.plot_pressure(20, 3)

while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    ns.step(ns.find_suitable_dt())
    
    if ch == 'p' or ch == 'P':
        plot_plots()
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)

kmkns.Keyboard.set_normal()

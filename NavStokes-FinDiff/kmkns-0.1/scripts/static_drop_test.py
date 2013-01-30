#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 's:l:w:m:', ['save=', 'load=', 'size=', 'method='])
loadfile = ""
savefile = ""
size = 40
method = "mdac"
for opt, val in opts:
    if opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val
    elif opt in ('-w', '--size'):
        size = int(val)
    elif opt in ('-m', '--method'):
        method = val

grid = Grid([-1.0, -1.0], [1.0, 1.0], [size, size])
ns = NavierStokes2D(grid, mu=0.1, mu_gas=1.6e-4,
    rho=1.0, rho_gas=1.3e-3, surface_tension_coeff=714.0, 
    gravity=[0.0, 0.0], curv_method=method)

ns.set_bc('west+east', u=0, v=0)
ns.set_bc('south+north', u=0, v=0)
ns.add_gas([-1.0, -1.0], [1.0, 1.0])
ns.add_drop([0.0, 0.0], 0.5)

level_contours = [0.3, 0.5, 0.7]
curv_contours = [-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.5, -2.3, -2.1, -2.05,
    -2.03, -2.01, -2.005, -2.003, -2.001, -2.0, -1.999, -1.997, -1.995, -1.99, -1.97,
    -1.95, -1.9, -1.7, -1.5, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]

kmkns.Keyboard.set_poll()

def plot_plots():
    ns.plot_level(level_contours, fig=1)
    #ns.plot_surface_tension(fig=2)
    ns.plot_velocity(fig=2)
    ns.plot_curvature(curv_contours, fig=3)
    ns.plot_pressure(20, fig=4)

    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "Max u-velocity: %g" % ns.u[:,1:-1].max()
    print "Max v-velocity: %g" % ns.v[1:-1,:].max()

plot_dt = 0.002

while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    dt = min(ns.find_suitable_dt(), plot_dt)
    ns.step(dt)
    
    if ch == 'p' or ch == 'P':
        plot_plots()                                
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)
        
kmkns.Keyboard.set_normal()

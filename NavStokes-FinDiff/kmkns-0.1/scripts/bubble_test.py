#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 's:l:w:', ['save=', 'load=', 'size='])
loadfile = ""
savefile = ""
size = 100
for opt, val in opts:
    if opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val
    elif opt in ('-w', '--size'):
        size = int(val)

grid = Grid([0.0, 0.0], [2.0, 4.0], [size, 2 * size])

ns = NavierStokes2D(grid, mu=0.002, mu_gas=3.2e-5,
    rho=1.0, rho_gas=0.0013, surface_tension_coeff=1.46, gravity=[0.0, -5.0])
    
ns.set_bc('west+east', u=0, v=0)
ns.set_bc('south+north', u=0, v=0)

ns.add_bubble([1.0, 1.0], 0.5)

kmkns.Keyboard.set_poll()

def plot_plots():
    ns.plot_level([0.5], fig=1)
    ns.plot_mass_centre_velocity(fig=2)
    ns.plot_gas_fraction(fig=3)
    #ns.plot_velocity(fig=4)
    print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
    print "# of cells > 0.5 = %d" % (ns.c > 0.5).sum()


while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    #dt = min(ns.find_suitable_dt(0.5), 0.02)
    #if (ns.t + dt >= next_plot_time) and (next_plot_time <= end_time):
    #    ns.step(dt)
    #    plot_plots()
    #    next_plot_time += plot_dt
    #else:
    #    ns.step(dt)
    
    ns.step(ns.find_suitable_dt())
        
    if ch == 'p' or ch == 'P':
        plot_plots()
    
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)
    
kmkns.Keyboard.set_normal()

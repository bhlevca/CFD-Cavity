#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 's:l:r:w:', ['save=','load=','Re=', 'size='])
loadfile = ""
savefile = ""
Re = 1.0
size = 25
for opt, val in opts:
    if opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val
    elif opt in ('-r', '--Re'):
        Re = float(val)
    elif opt in ('-w', '--size'):
        size = float(val)

grid = Grid([0.0, 0.0], [3.0, 1.0], [3 * size, 2 * size + 1])
ns = NavierStokes2D(grid, nu=1.0/Re)

ns.set_bc('north+south', u=0, v=0)
ns.set_bc('west', u=1, v=0)
ns.set_bc('east', dudn=0, dvdn=0)

kmkns.Keyboard.set_poll()
while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break
        
    ns.step(ns.find_suitable_dt(), beta=0.31337, gamma=0.0)

    if ch == 'p' or ch == 'P':
        ns.plot_velocity(1)
        ns.plot_pressure(20, 2)
        #ns.plot_velocity_profile(0, [0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0], 3)
        #ns.plot_velocity_component(0, None, 4)
        #ns.plot_divergence(None, 5)
        #ns.plot_vorticity(None, 6)
        #ns.plot_stream(20, 7)

        print "Time: %f" % ns.t
        print "Iteration #%i" % ns.it
        mid = ns.u.shape[1]/2
        print "Centre velocities (u):", ns.u[-1,mid-2:mid+3]
        print "Max velocity (u):", ns.u[-1,mid]
        print "Error (u):", 1.5 - ns.u[-1,mid]
            
    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)
    
kmkns.Keyboard.set_normal()

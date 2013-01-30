#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:r:s:l:', ['size=', 'Re=', 'save=', 'load='])
size = 10
Re = 100.0
savefile = ""
loadfile = ""
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)
    elif opt in ('-r', '--Re'):
        Re = float(val)
    elif opt in ('-s', '--save'):
        savefile = val
    elif opt in ('-l', '--load'):
        loadfile = val

grid = Grid([-2.0, -2.0], [20.0, 2.1], [22 * size, 41 * size / 10])

kmkns.Keyboard.set_poll()

ns = NavierStokes2D(grid, mu=1.0/Re)
ns.set_bc('south+north', u=0, v=0)
inlet_vel = lambda x, y, t: 6.0 * (y + 2.0) * (2.1 - y) / 4.1**2
ns.set_bc('west', u=inlet_vel, v=0)
ns.set_bc('east', dudn=0, dvdn=0)
ns.set_obstacles(lambda x, y: (x**2 + y**2) < 0.25)
ns.add_tracer([-3.0, -2.0], [-1.0, 0.0])

def plot_plots():
    print "t=%g" % ns.t
    ns.plot_stream(20, fig=1)
    ns.plot_velocity(fig=2)
    ns.plot_tracer(fig=3)

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

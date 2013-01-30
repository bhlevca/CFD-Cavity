#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
savefile = ""
loadfile = ""
size = 50
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = float(val)

kmkns.Keyboard.set_poll()

grid = Grid([-1.0, -1.0], [1.0, 1.0], [size, size])
ns = NavierStokes2D(grid, mu_liquid=1.0, mu_gas=1.0, rho_liquid=1.0, rho_gas=1.0)

u_bc = lambda x, y, t: -y
v_bc = lambda x, y, t: x
u_ic = lambda x, y: -y
v_ic = lambda x, y: x
ns.set_ic(u=u_ic, v=v_ic)
ns.set_bc('north+south+west+east', u=u_bc, v=v_bc)

ns.add_bubble([0, 0.5], 0.3)
ns.add_liquid([-0.05, -1.0], [0.05, 0.7])

def show_info():
    print "Iteration #%i" % ns.it
    print "Time: %f" % ns.t
    ns.plot_level([0.5], 1)

while ns.t < 2.0 * pi:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break

    ns.step(ns.find_suitable_dt())

    if ch == 'p' or ch == 'P':
        show_info()

    if (ch == 'l' or ch == 'L') and loadfile:
        ns.load(loadfile)
    if (ch == 's' or ch == 'S') and savefile:
        ns.save(savefile)

kmkns.Keyboard.set_normal()

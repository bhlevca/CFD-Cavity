#!/usr/bin/env python
from kmkns.NavierStokes import *
import kmkns.Keyboard

grid = Grid([0.0, 0.0], [1.0, 1.0], [32, 32])
ns = NavierStokes2D(grid, mu_liquid=1.0, mu_gas=1e-2,
    rho_liquid=1000.0, rho_gas=1.0, surface_tension_coeff=0.0, gravity=[0.0, -100.0])

ns.set_bc('north+south+west+east', u=0, v=0)
ns.add_gas([0.0, 0.0], [1.0, 1.0])
ns.add_drop([0.5, 0.5], 0.25)

print "Press 'q' to quit and 'p' to plot."

kmkns.Keyboard.set_poll()
while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break        
    ns.step(ns.find_suitable_dt())
    if ch == 'p' or ch == 'P':
        ns.plot_level([0.5])
        print 't=%g' % ns.t
    
kmkns.Keyboard.set_normal()

print 'Bye, bye!'

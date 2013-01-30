#!/usr/bin/env python
from kmkns.NavierStokes import *
import kmkns.Keyboard

grid = Grid([0.0, 0.0], [3.0, 1.0], [90, 30])
ns = NavierStokes2D(grid, mu=1.0, rho=1.0)

ns.set_bc('north+south', u=0, v=0)
ns.set_bc('west', u=1, v=0)
ns.set_bc('east', dudn=0, dvdn=0)

print "Press 'q' to quit and 'p' to plot."

kmkns.Keyboard.set_poll()
while True:
    ch = kmkns.Keyboard.get_key()
    if ch == 'q' or ch == 'Q':
        break        
    ns.step(ns.find_suitable_dt())
    if ch == 'p' or ch == 'P':
        ns.plot_velocity(1)
    
kmkns.Keyboard.set_normal()

print 'Bye, bye!'

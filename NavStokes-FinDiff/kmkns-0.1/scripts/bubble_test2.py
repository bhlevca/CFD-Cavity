#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

filename = "bubble_test2"
sizes = [25, 50, 100, 200]
end_time = 0.5

kmkns.Keyboard.set_poll()

for size in sizes:
    grid = Grid([0.0, 0.0], [2.0, 4.0], [size, 2 * size])
    ns = NavierStokes2D(grid, mu=0.002, mu_gas=3.2e-5,
        rho=1.0, rho_gas=0.0013, surface_tension_coeff=1.46, gravity=[0.0, -5.0])

    ns.set_bc('west+east+south+north', u=0, v=0)
    ns.add_bubble([1.0, 1.0], 0.5)

    ns.load("%s_%i.dat" % (filename, size))

    def plot_plots(save):
        ns.plot_level([0.5], fig=1)
        if save:
            if viz.backend == 'gnuplot':
                g = viz.get_backend()
                g("set size ratio 2")
                g("set xrange [0:2]")
                g("set yrange [0:4]")
                g("set rmargin 10")
                g("set lmargin 10")
                g("set bmargin 10")
                g("set tmargin 10")
                g("unset key")
                g("replot")
                g.hardcopy("bubble_%il.eps" % size)
            else:
                viz.hardcopy("bubble_%il.eps" % size)

        ns.plot_gas_centre_velocity(fig=2)
        if save:
            if viz.backend == 'gnuplot':
                g = viz.get_backend()
                g("unset key")
                g("replot")
                g.hardcopy("bubble_%iv.eps" % size)
            else:
                viz.hardcopy("bubble_%iv.eps" % size)

        ns.plot_gas_fraction(fig=3)
        viz.legend("Area of gas/Total area", "Number of gas cells/Total number of cells")
        if save:
            if viz.backend == 'gnuplot':
                g = viz.get_backend()
                g("replot")
                g.hardcopy("bubble_%im.eps" % size)
            else:
                viz.hardcopy("bubble_%im.eps" % size)

        #ns.plot_mass_centre(fig=4)
        #ns.plot_velocity(fig=4)
        print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
        print "# of cells > 0.5 = %d" % (ns.c[1:-1,1:-1] > 0.5).sum()
        print "Sum of colour function = %g" % ns.c[1:-1,1:-1].sum()

    while ns.t < end_time:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break

        dt = ns.find_suitable_dt()
        if (ns.t + dt >= end_time):
            ns.step(end_time - ns.t)
        else:
            ns.step(dt)

        if ch == 'p' or ch == 'P':
            plot_plots(False)

        if (ch == 'l' or ch == 'L'):
            ns.load("%s_%i.dat" % (filename, size))
        if (ch == 's' or ch == 'S'):
            ns.save("%s_%i.dat" % (filename, size))

    plot_plots(True)
    ns.save("%s_%i.dat" % (filename, size))

dummy_a = 1
dummy_b = 1
for dummy_loop in xrange(200000):
    dummy_c = dummy_a + dummy_b
    dummy_a = dummy_b
    dummy_b = dummy_c

kmkns.Keyboard.set_normal()

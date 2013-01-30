#!/usr/bin/env python

import sys, getopt
from kmkns.NavierStokes import *
from numpy import *
import kmkns.Keyboard

opts, args = getopt.getopt(sys.argv[1:], 'w:', ['size='])
size = 96
surf_coeff = 0.357
for opt, val in opts:
    if opt in ('-w', '--size'):
        size = int(val)

level_contours = [0.5]
curv_contours = [-5.0, -4.3, -4.1, -4.03, -4.01, -4.003, -4.001, 
    -4.0, -3.999, -3.997, -3.99, -3.97, -3.9, -3.7, -3.0]
radius = 0.25

grid = Grid([0.0, 0.0], [1.0, 1.0], [size, size])

kmkns.Keyboard.set_poll()

for method in ("mdac", "dac", "csf"):
    ns = NavierStokes2D(grid, mu=1.0, mu_gas=1.0,
        rho=4.0, rho_gas=4.0, surface_tension_coeff=surf_coeff, 
        gravity=[0.0, 0.0], curv_method=method)
    ns.set_bc('west+east', u=0, v=0)
    ns.set_bc('south+north', u=0, v=0)
    ns.add_gas([0.0, 0.0], [1.0, 1.0])
    ns.add_drop([0.5, 0.5], radius)

    filename = "static_%s.dat" % method
    ns.load(filename)

    def plot_plots():
        ns.plot_level(level_contours, fig=1)
        #ns.plot_surface_tension(fig=2)
        ns.plot_velocity(fig=2)
        ns.plot_curvature(curv_contours, fig=3)
        ns.plot_pressure(20, fig=4)

        print "Plotted at iteration #%d, time=%f" % (ns.it, ns.t)
        print "Max u-velocity: %g" % ns.u[:,1:-1].max()
        print "Max v-velocity: %g" % ns.v[1:-1,:].max()
        pressure = ns.p[1 + size / 2, 1 + size / 2]
        error = abs(surf_coeff / radius - pressure)
        print "Centre pressure: %g" % pressure
        print "Expected pressure: %g" % (surf_coeff / radius)
        print "Pressure error: %g" % error
        print "Relative error: %g" % (error / abs(surf_coeff / radius))

    dt = 1e-5
    plot_dt = 1e-4
    next_plot_time = 1e-4
    end_time = 200.5e-5
    while next_plot_time <= ns.t:
        next_plot_time += plot_dt

    while ns.t < end_time:
        ch = kmkns.Keyboard.get_key()
        if ch == 'q' or ch == 'Q':
            break

        if (ns.t + dt >= next_plot_time) and (next_plot_time <= end_time):
            ns.step(next_plot_time - ns.t)
            plot_plots()

            next_plot_time += plot_dt
            if next_plot_time >= end_time:
                # save plots
                def save_plot(filename, shrink):
                    if viz.backend == 'gnuplot':
                        g = viz.get_backend()
                        g("set size square")
                        g("set xrange [0:1]")
                        g("set yrange [0:1]")
                        if shrink:
                            g("set rmargin 12")
                            g("set lmargin 12")
                            g("set bmargin 12")
                            g("set tmargin 12")
                        g("replot")
                        g.hardcopy(filename)
                    else:
                        viz.hardcopy(filename)

                ns.plot_level(level_contours, fig=1)
                save_plot("static_level_%s.eps" % method, True)
                ns.plot_velocity(fig=2)
                save_plot("static_velocity_%s.eps" % method, False)
                ns.plot_curvature(curv_contours, fig=3)
                save_plot("static_curvature_%s.eps" % method, True)
                ns.plot_pressure(20, fig=4)
                save_plot("static_pressure_%s.eps" % method, True)
        else:
            ns.step(dt)

        if ch == 'p' or ch == 'P':
            plot_plots()
        if ch == 'l' or ch == 'L':
            ns.load(filename)
        if ch == 's' or ch == 'S':
            ns.save(filename)

    if ns.t < end_time:
        break
        
    #ns.save(filename)

    dummy_a = 1
    dummy_b = 1
    for dummy_loop in xrange(200000):
        dummy_c = dummy_a + dummy_b
        dummy_a = dummy_b
        dummy_b = dummy_c

kmkns.Keyboard.set_normal()

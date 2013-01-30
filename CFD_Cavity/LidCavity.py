'''
Created on Apr 7, 2012

@author: bogdan
'''

import numpy
import SIMPLE

#import the graphi library
import matplotlib.pyplot as plt
from streamplot import streamplot

class LidCavity(object):
    '''
    classdocs
    '''


    def __init__(self, scheme, debug):
        '''
        Constructor
        '''
        self.W = 0
        self.E = 1
        self.N = 2
        self.S = 3
        self.P = 4


        self.CPP = True
        #domain
        self.Lx = 1.0
        self.Ly = 1.0
        self.Nx = 25  #25  -for QUICK #129
        self.Ny = 25  #25  -for QUICK #i29

        #fluxes for u, v, p
        self.u0 = 1.0
        self.FiU = [0, 0, self.u0, 0 ]
        self.FiV = [0, 0, 0, 0 ]
        #self.FiP = [0, 0, 0, 0 ]

        self.ue_wall = numpy.zeros(self.Nx) # horizontal wall velocity
        self.uw_wall = numpy.zeros(self.Nx) # horizontal wall velocity
        self.vn_wall = numpy.zeros(self.Ny) # vertical wall velocity
        self.vs_wall = numpy.zeros(self.Ny) # vertical wall velocity

        for j in range(0, self.Ny):
            self.ue_wall[j] = self.FiU[self.E]
            self.uw_wall[j] = self.FiU[self.W]

        for i in range(0, self.Nx):
            self.vn_wall[i] = self.FiV[self.N]
            self.vs_wall[i] = self.FiV[self.S]


        #UD  GOOD run 25x25, 50x50 , 129x129 urf_p =0.1 urf_uv=0.6 all corrections in m corrections 41 iter. Re=100 rho = 10
        #QUICK GOOD run 25x25, urf_p =0.06 urf_uv=0.8 all corrections in m corrections 41 iter. Re=100 rho = 10


        #initial conditions at boundary
        self.Re = 100
        self.rho = 10
        self.miuX = self.rho * self.FiU[self.N] * self.Lx / self.Re
        self.miuY = self.rho * self.FiU[self.N] * self.Ly / self.Re
        self.simple = 0
        self.scheme = scheme
        self.debug = debug






    def solve(self):
        self.simple = SIMPLE.SIMPLE(self, self.scheme, self.debug)
        u, v = self.simple.solve()
        return [u, v]

    def plotStreams(self, u, v):
        x = numpy.linspace(0, self.Lx , self.Nx)
        y = numpy.linspace(0, self.Ly, self.Ny)
        speed = numpy.sqrt(u * u + v * v)
        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')
        plt.subplot(111)
        #streamplot(x, y, u, v, density = (1, 1), INTEGRATOR = 'RK4', color = u, linewidth = 1 * speed / speed.max())
        streamplot(x, y, u, v, density = (1, 1), INTEGRATOR = 'RK4', color = 'b', linewidth = 10 * speed / speed.max())

    def plotHeatMap(self, u, v):
        x = numpy.linspace(0, self.Lx , self.Nx)
        y = numpy.linspace(0, self.Ly, self.Ny)
        speed = numpy.sqrt(u * u + v * v)
        X, Y = numpy.meshgrid(x, y)                             #create the grid
        #plot a heatmap
        fig2 = plt.figure(2, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
        ax = fig2.add_subplot(111)
        im = ax.pcolor(X, Y, speed)                             # draw the heatmap
        fig2.colorbar(im)                                        # add the legend on a colour bar

        #superimpose contours
        cs = ax.contour(X, Y, speed)                          # draw the contour 
        plt.clabel(cs, inline = 1, fontsize = 10)               # draw the contour

        #draw only the contours
        fig3 = plt.figure(3, facecolor = 'w', edgecolor = 'k')
        ax2 = fig3.add_subplot(111)
        cs = ax2.contour(X, Y, speed)                          # draw the contour
        plt.clabel(cs, inline = 1, fontsize = 10)                # create the label on the contour



    def plotCentreLinesProfiles(self, u, v):
        x = numpy.linspace(0, self.Lx , self.Nx)
        y = numpy.linspace(0, self.Ly, self.Ny)
        #take center values
        uc = u[:, self.Ny / 2]
        vc = v [self.Nx / 2, :]
        fig4 = plt.figure(4, facecolor = 'w', edgecolor = 'k')
        ax4 = fig4.add_subplot(121)
        ax4.plot(y, uc)                          # draw the contour
        ax4.set_ylabel("u- gravity centre")
        ax4.set_xlabel("length")
        ax4.grid(b = True)

        ax5 = fig4.add_subplot(122)
        ax5.plot(y, vc)                          # draw the contour
        ax5.set_ylabel("v - gravity centre")
        ax5.set_xlabel("length")
        ax5.grid(b = True)

    def showGraph(self):
        plt.show()



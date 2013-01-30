'''
Created on Feb 6, 2012

@author: bogdan
'''

import numpy
from thomas import *
import matplotlib.pyplot as plt






class FinVol_2D_Diffusion_Ex7_2:

    def __init__(self, Lx, Ly, Gx, Gy, FiN, Nx, Ny, qw, d, debug = False):
        self.debug = debug
        self.W = 0
        self.E = 1
        self.N = 2
        self.S = 3

        self.Gamma = numpy.zeros(4)
        self.A = numpy.zeros(4)

        self.deltaX = None
        self.deltaY = None

        #iterative error accepted
        self.EPSILON = 1e-4

        self.FiN = FiN
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        self.QW = qw
        self.d = d

        self.initNodesEx()
        self.calculateGammaEx(Gx, Gy)
        self.calculateDeltaEx()
        self.calculateAEx()

    def calculateGammaEx(self, GammaX, GammaY):
        self.Gamma[self.W] = self.Gamma[self.E] = GammaY
        self.Gamma[self.N] = self.Gamma[self.S] = GammaX

    def calculateAEx(self):
        self.A[self.W] = self.A[self.E] = self.deltaY * self.d
        self.A[self.N] = self.A[self.S] = self.deltaX * self.d

    def calculateDeltaEx(self):
        self.deltaX = self.Lx / self.Nx
        self.deltaY = self.Ly / self.Ny



    def initNodesEx(self):
        '''
        '''
        #init the nodes
        self.aW = numpy.zeros((self.Ny, self.Nx))
        self.aE = numpy.zeros((self.Ny, self.Nx))
        self.aN = numpy.zeros((self.Ny, self.Nx))
        self.aS = numpy.zeros((self.Ny, self.Nx))
        self.aP = numpy.zeros((self.Ny, self.Nx))

        #init tself.SpWhe sources
        self.SuN = numpy.zeros((self.Ny, self.Nx))
        self.SuS = numpy.zeros((self.Ny, self.Nx))
        self.SuE = numpy.zeros((self.Ny, self.Nx))
        self.SuW = numpy.zeros((self.Ny, self.Nx))

        self.SpN = numpy.zeros((self.Ny, self.Nx))
        self.SpS = numpy.zeros((self.Ny, self.Nx))
        self.SpE = numpy.zeros((self.Ny, self.Nx))
        self.SpW = numpy.zeros((self.Ny, self.Nx))

        #init with zero the solution matrix     
        self.Fi = numpy.zeros((self.Ny, self.Nx))
        self.FiOld = numpy.zeros((self.Ny, self.Nx))

        #the TDMA coeficients. (We calculate vertical lines (columns = > allocate the number of horizontal lines)
        self.alp = numpy.zeros(self.Ny)
        self.bet = numpy.zeros(self.Ny)
        self.D = numpy.zeros(self.Ny)

        self.C = numpy.zeros(self.Ny)

    #end initNodes


    def calculateSourcesEx(self):
        '''
        '''
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #west
                if i == 0:
                    self.SuW[j, i] = self.QW * self.A[self.W]
                    self.SpW[j, i] = 0
                else:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0

                #east
                if i == self.Nx - 1:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0
                else:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0

                #north has a flux h(Fi-Fi0)
                if j == self.Ny - 1:
                    self.SuN[j, i] = 2 * self.Gamma[self.N] * self.A[self.N] / self.deltaY * self.FiN
                    self.SpN[j, i] = -2 * self.Gamma[self.N] * self.A[self.N] / self.deltaY
                else:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0
                #south
                if j == 0:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
        #end j
        if self.debug == True:
            print "Su:", self.SuE + self.SuW + self.SuS + self.SuN
            print "Sp:", self.SpE + self.SpW + self.SpS + self.SpN

    #end calculateSources 

    def calculateCoefficientsEx(self):
        '''
        '''
        self.calculateSourcesEx()

        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                #west
                if i == 0:
                    self.aW[j, i] = 0
                else:
                    self.aW[j, i] = self.Gamma[self.W] / self.deltaX * self.A[self.W]

                #east
                if i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    self.aE[j, i] = self.Gamma[self.E] / self.deltaX * self.A[self.E]

                #north has a flux h(Fi-Fi0)
                if j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    self.aN[j, i] = self.Gamma[self.N] / self.deltaY * self.A[self.N]
                #south
                if j == 0:
                    self.aS[j, i] = 0
                else:
                    self.aS[j, i] = self.Gamma[self.S] / self.deltaY * self.A[self.S]

                #ap coefficient
                Sp = self.SpW[j, i] + self.SpE[j, i] + self.SpS[j, i] + self.SpN[j, i]
                self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] - Sp

        #END for i


    #end calculateCoefficients





    def calculateTDMACoefficientsEx(self, i):
        '''
        book pag 220
        Apply TDMA S to N sweeping W to E
        The discretisation equation is given by
       
        In the book they have it reversed "j" is for lines and "i" for columns 
        '''
        #calculate on each vertical  from S -> N 
        for j in range(0, self.Ny):

            #set boundary condition for E and W of line i if i=0 or i=N
            self.alp[j] = self.aN[j, i].copy()
            self.bet[j] = self.aS[j, i].copy()
            self.D[j] = self.aP[j, i].copy()

            #the free term
            Su = self.SuW[j, i] + self.SuE[j, i] + self.SuS[j, i] + self.SuN[j, i]
            if i == 0 :
                FiW = 0
            else:
                FiW = self.Fi[j, i - 1]

            if i == self.Nx - 1:
                FiE = 0
            else :
                FiE = self.Fi[j, i + 1]

            self.C[j] = self.aW[j, i] * FiW + self.aE[j, i] * FiE + Su


        #end for j


    #end calculateTDMACoefficients

    def callTDMAEx(self):
        #solve with TDMA for this i vertical line S-N
        n = self.D.size
        maxiter = 500
        iter = 0

        while  maxiter > iter :
            #Swipe  from W to E
            for i in range(0, self.Nx):

                #calculate the TDMA coefficients for column i
                self.calculateTDMACoefficientsEx(i)

                if self.debug == True:
                    print "beta:", self.bet
                    print "D", self.D
                    print "alp", self.alp
                    print "C", self.C


            #for i in range(1, self.Nx - 1):

                #Provide the solution vector for the C++ 
                #x = self.Fi[:, i]
                #d = cls.solveTDMA(x, n)
                x = thomas(n, -self.bet[1:], self.D, -self.alp[:-1], self.C)
                self.Fi[:, i] = x.copy()
            #end i
            #test accuracy and exit condition 
            flat = self.Fi[1] - self.FiOld[1]
            dx = math.sqrt(numpy.dot(flat, flat))

            print "iteration # %d, dx=%f" % (iter, dx)
            #print "Fi:", self.Fi

            if dx < self.EPSILON :
                return
            iter += 1
            #copy the old values
            self.FiOld = self.Fi.copy()

            if maxiter < iter:
                print "Max iterations exceeded => did not converge"
                return
        #end while 
    #end callTDMA

    def plotFiEx(self):

        #create the figure and axes object
        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')
        ax = fig.add_subplot(111)

        # Creating the grid of coordinates x,y 
        x = numpy.arange(0., self.Lx, self.deltaX)
        y = numpy.arange(0., self.Ly, self.deltaY)
        X, Y = numpy.meshgrid(x, y)

        #ax.contour(X, Y, self.Fi)
        im = ax.pcolor(X, Y, self.Fi)
        fig.colorbar(im)
        plt.show()

    def solveEx(self):
        self.calculateCoefficientsEx()

        if self.debug == True:
            print "aW coefficients:"
            print self.aW
            print "aE coefficients:"
            print self.aE
            print "aS coefficients:"
            print self.aS
            print "aN coefficients:"
            print self.aN
            print "aP coefficients:"
            print self.aP

        self.callTDMAEx()
        print self.Fi
        self.plotFiEx()

    #end solve

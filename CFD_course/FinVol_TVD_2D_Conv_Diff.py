'''
Created on Feb 6, 2012

@author: bogdan hlevca 995151213
'''

#import numerical Python for Matrix operations
import numpy

#import the TDMA module
from thomas import *

#import the graphi library
import matplotlib.pyplot as plt


#Class definition
class FinVol_TVD_2D_Conv_Diff:
    '''
    Solve the equation:
    d/dx(ro*u*F) + d/dy(ro*v*F) = d/dx(Gamma* dF/dx)+d/dy(Gamma* dF/dy) 
    
    where 
        d is the partial difference d_rond
        F - stands for the property Fi
    '''

    #constructor object initializes the class arrays and calculates preliminary coefficients
    def __init__(self, Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, scheme, debug = False):
        '''
        Constructor
        '''

        self.debug = debug
        self.W = 0
        self.E = 1
        self.N = 2
        self.S = 3

        self.rep = 0 #re+
        self.rem = 1 #re-
        self.rwp = 2 #rw+
        self.rwm = 3 #rw-
        self.rnp = 4 #rn+
        self.rnm = 5 #rn-
        self.rsp = 6 #rs+
        self.rsm = 7 #rs-

        self.Gamma = numpy.zeros(4)
        self.A = numpy.zeros(4)
        self.F = numpy.zeros(4)
        self.Df = numpy.zeros(4)

        self.deltaX = None
        self.deltaY = None

        #iterative error accepted
        self.EPSILON = 1e-2
        self.maxiter = 2000

        self.Fi0E = FiE
        self.Fi0N = FiN
        self.Fi0W = FiW
        self.Fi0S = FiS

        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        self.scheme = scheme
        self.u = u
        self.v = v
        self.rho = rho



        self.initNodes()               #initialize node coefficients
        self.calculateGamma(Gx, Gy)    #calculate preliminary information, Gamma, deltaX and Y and the control volume area 
        self.calculateDelta()
        self.calculateA()
        self.calculate_F()
        self.calculate_Df()

    #end __init__

    def initNodes(self):
        '''
        initialize node coefficients and TDMA coefficients
        '''
        #init the nodes
        self.aW = numpy.zeros((self.Ny, self.Nx))
        self.aE = numpy.zeros((self.Ny, self.Nx))
        self.aN = numpy.zeros((self.Ny, self.Nx))
        self.aS = numpy.zeros((self.Ny, self.Nx))
        self.aP = numpy.zeros((self.Ny, self.Nx))

        #init Su sources
        self.SuN = numpy.zeros((self.Ny, self.Nx))
        self.SuS = numpy.zeros((self.Ny, self.Nx))
        self.SuE = numpy.zeros((self.Ny, self.Nx))
        self.SuW = numpy.zeros((self.Ny, self.Nx))

        #init Sp sources
        self.SpN = numpy.zeros((self.Ny, self.Nx))
        self.SpS = numpy.zeros((self.Ny, self.Nx))
        self.SpE = numpy.zeros((self.Ny, self.Nx))
        self.SpW = numpy.zeros((self.Ny, self.Nx))

        #Deferred correction
        self.SuDC = numpy.zeros((self.Ny, self.Nx))

        #init with zero the solution matrix     
        self.Fi = numpy.zeros((self.Ny, self.Nx))
        self.FiOld = numpy.zeros((self.Ny, self.Nx))



        #the TDMA coeficients. (We calculate vertical lines (acolumns = > allocate the number of horizontal lines)
        self.alp = numpy.zeros(self.Ny)
        self.bet = numpy.zeros(self.Ny)
        self.D = numpy.zeros(self.Ny)
        self.C = numpy.zeros(self.Ny)

    def calculateGamma(self, GammaX, GammaY):
        '''
        calculate preliminary information, Gamma
        '''
        self.Gamma[self.W] = self.Gamma[self.E] = GammaY
        self.Gamma[self.N] = self.Gamma[self.S] = GammaX

    def calculateA(self):
        '''
        calculate preliminary information, control volume area
        '''
        self.A[self.W] = self.A[self.E] = self.deltaY
        self.A[self.N] = self.A[self.S] = self.deltaX

    def calculateDelta(self):
        '''
        calculate preliminary information, deltaX and Y 
        '''
        self.deltaX = self.Lx / self.Nx
        self.deltaY = self.Ly / self.Ny

    def calculate_F(self):
       '''steady flow F is constant
       
       '''
       self.F[self.W] = self.rho * self.u * self.A[self.W]
       self.F[self.E] = self.rho * self.u * self.A[self.E]
       self.F[self.S] = self.rho * self.v * self.A[self.S]
       self.F[self.N] = self.rho * self.v * self.A[self.N]
       self.calculate_alpha()

    #end Calculate_F

    def calculate_Df(self):
        '''
        '''
        self.Df[self.W] = self.Gamma[self.W] * self.A[self.W] / self.deltaX
        self.Df[self.E] = self.Gamma[self.E] * self.A[self.E] / self.deltaX
        self.Df[self.S] = self.Gamma[self.S] * self.A[self.S] / self.deltaY
        self.Df[self.N] = self.Gamma[self.N] * self.A[self.N] / self.deltaY
    #end calculate_D

    def calculate_alpha(self):
        if self.F[self.W] >= 0:
            self.Alpha_w = 1
        else :
            self.Alpha_w = 0

        if self.F[self.E] >= 0:
            self.Alpha_e = 1
        else :
            self.Alpha_e = 0

        if self.F[self.S] >= 0:
            self.Alpha_s = 1
        else :
            self.Alpha_s = 0

        if self.F[self.N] >= 0:
            self.Alpha_n = 1
        else :
            self.Alpha_n = 0
    #end calculate_alpha






    def calculateSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions
        
        For the Assignment 3 all the sources  = 0 since there is no flux on the borders
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #west
                if i == 0:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0
                elif i == 1 :
                    self.SuW[j, i] = 0
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

                #south
                if j == 0:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
                elif j == 1:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0

                #north
                if j == self.Ny - 1:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0
                else:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0

        #end j
        if self.debug == True:
                print "Su:", self.SuE + self.SuW + self.SuS + self.SuN
                print "Sp:", self.SpE + self.SpW + self.SpS + self.SpN

    #end calculateSources 

    def setBoundaryConditions(self):
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):


                #north
                if j == self.Ny - 1:
                    self.Fi[j, i] = self.Fi0N

                #west
                if i == 0:
                    self.Fi[j, i] = self.Fi0W

                #east
                if i == self.Nx - 1:
                    self.Fi[j, i] = self.Fi0E

                #south
                if j == 0:
                    self.Fi[j, i] = self.Fi0S

    #end set BoundaryConditions


    def PSI(self, r):
        psi = 0
        if self.scheme == 'UD':
            psi = 0
        elif self.scheme == 'CD':
            psi = 1
        elif self.scheme == 'LUD':
            psi = r
        elif self.scheme == 'QUICK':
            psi = (3 + r) / 4
        elif self.scheme == 'VanLeer':
            psi = (r + abs(r)) / (1 + r)

        return psi
    #end PSI

    def FiW(self, j, i):
        #FiW
        if i == 0:
            FiW = self.Fi0W
        else:
            FiW = self.Fi[j, i - 1]
        return FiW

    def FiE(self, j, i):
        if i == self.Nx - 1:
            FiE = self.Fi0E
        else:
            FiE = self.Fi[j, i + 1]
        return FiE

    def FiS(self, j, i):
        if j == 0:
            FiS = self.Fi0S
        else:
            FiS = self.Fi[j - 1, i]
        return FiS

    def FiN(self, j, i):
        if j == self.Ny - 1:
            FiN = self.Fi0N
        else:
            FiN = self.Fi[j + 1, i]
        return FiN

    def FiWW(self, j, i):
        if i == 0:
            FiWW = 2 * self.Fi0W - self.Fi[j, i]  #4 * self.Fi0W - 3 * self.Fi[j, i]
        elif i == 1:
            FiWW = self.Fi0W #2 * self.Fi0W - self.FiW(j, i)
        else:
            FiWW = self.Fi[j, i - 2]
        return FiWW

    def FiEE(self, j, i):
        if i == self.Nx - 1  :
            FiEE = 2 * self.Fi0E - self.Fi[j, i]#4 * self.Fi0E - 3 * self.Fi[j, i]
        elif  i == self.Nx - 2:
            FiEE = self.Fi0E
        else:
            FiEE = self.Fi[j, i + 2]
        return FiEE

    def FiSS(self, j, i):
        if j == 0:
            FiSS = 2 * self.Fi0S - self.Fi[j, i]#4 * self.Fi0S - 3 * s2 * self.Fi[j, i]
        elif j == 1:
            FiSS = self.Fi0S
        else:
            FiSS = self.Fi[j - 2, i]
        return FiSS

    def FiNN(self, j, i):
        if j == self.Ny - 1  :
            FiNN = 2 * self.Fi0N - self.Fi[j, i]#4 * self.Fi0N - 3 * self.Fi[j, i]
        elif j == self.Ny - 2:
            FiNN = self.Fi0N#2 * self.Fi0N - self.FiN(j, i)
        else:
            FiNN = self.Fi[j + 2, i]
        return FiNN




    def r(self, t, j, i):
        '''
            i - columns
            j - lines
            rep == re+
            rem == re-
            rwp == rw+
            rwm == rw-
        '''
        rr = "unproc"
        #if u >= 0:
        if t == self.rep:
            #rep = (FiP - FiW) / (FiE - FiP)
            if (self.FiE(j, i) - self.Fi[j, i]) == 0:
                return 0.0
            else:
                rr = (self.Fi[j, i] - self.FiW(j, i)) / (self.FiE(j, i) - self.Fi[j, i]) #verified
        elif t == self.rwp:
            #rwp = (FiW - FiWW) / (FiP - FiW)
            if (self.Fi[j, i] - self.FiW(j, i)) == 0:
                return 0.0
            else:
                rr = (self.FiW(j, i) - self.FiWW(j, i)) / (self.Fi[j, i] - self.FiW(j, i))
        # if u < 0 
        elif t == self.rem:
            #rem = (FiEE - FiE) / (FiE - FiP)
            if (self.FiE(j, i) - self.Fi[j, i]) == 0:
                return 0.0
            else:
                rr = (self.FiEE(j, i) - self.FiE(j, i)) / (self.FiE(j, i) - self.Fi[j, i]) #verified
        elif t == self.rwm:
            #rwm = (FiE - FiP) / (FiP - FiW)
            if (self.Fi[j, i] - self.FiW(j, i)) == 0:
                return 0.0
            else:
                rr = (self.FiE(j, i) - self.Fi[j, i]) / (self.Fi[j, i] - self.FiW(j, i))

        # S -> N
        #if u >= 0:
        if t == self.rnp:
            #rnp = (FiP - FiS) / (FiN - FiP)
            if (self.FiN(j, i) - self.Fi[j, i]) == 0:
                return 0.0
            else:
                 rr = (self.Fi[j, i] - self.FiS(j, i)) / (self.FiN(j, i) - self.Fi[j, i]) #verified
        elif t == self.rsp:
            #rsp = (FiS - FiSS) / (FiP - FiS)
            if (self.Fi[j, i] - self.FiS(j, i)) == 0:
                return 0.0
            else:
                rr = (self.FiS(j, i) - self.FiSS(j, i)) / (self.Fi[j, i] - self.FiS(j, i))
        # if u < 0 
        elif t == self.rnm:
            #rnm = (FiNN - FiN) / (FiN - FiP)
            if (self.FiN(j, i) - self.Fi[j, i]) == 0:
                return 0.0
            else:
                rr = (self.FiNN(j, i) - self.FiN(j, i)) / (self.FiN(j, i) - self.Fi[j, i]) #verified
        elif t == self.rsm:
            #rsm = (FiN - FiP) / (FiP - FiS)
            if (self.Fi[j, i] - self.FiS(j, i)) == 0:
                return 0.0
            else:
                rr = (self.FiN(j, i) - self.Fi[j, i]) / (self.Fi[j, i] - self.FiS(j, i))

        if rr == "unproc":
            print "Wrong 'r' value passed. Exiting"
            raise Exception("Wrong 'r' value passed. Exiting !!")
        return rr
    #end r



    def calculateTVDcoefficients(self):
        '''
            2D "a" coeficients for TVD are at page 175 in the Versteeg book 
            
        '''
        self.calculateSources()


        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                if i == 0:
                    self.aW[j, i] = 0.0
                else:
                    self.aW[j, i] = self.Df[self.W] + max(self.F[self.W], 0)

                if i == self.Nx - 1:
                    self.aE[j, i] = 0.0
                else:
                    self.aE[j, i] = self.Df[self.E] + max(-self.F[self.E], 0)

                if j == 0:
                    self.aS[j, i] = 0.0
                else:
                    self.aS[j, i] = self.Df[self.S] + max(self.F[self.S], 0)

                if j == self.Ny - 1:
                    self.aN[j, i] = 0.0
                else:
                    self.aN[j, i] = self.Df[self.N] + max(-self.F[self.N], 0)

                Sp = self.SpW[j, i] + self.SpE[j, i] + self.SpS[j, i] + self.SpN[j, i]

                self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp

                #TVD deferred correction source term

                et = 0.5 * self.F[self.E] * ((1 - self.Alpha_e) * self.PSI(self.r(self.rem, j, i)) - self.Alpha_e * self.PSI(self.r(self.rep, j, i))) * (self.FiE(j, i) - self.Fi[j, i])
                wt = 0.5 * self.F[self.W] * (self.Alpha_w * self.PSI(self.r(self.rwp, j , i)) - (1 - self.Alpha_w) * self.PSI(self.r(self.rwm, j, i))) * (self.Fi[j, i] - self.FiW(j, i))
                nt = 0.5 * self.F[self.N] * ((1 - self.Alpha_n) * self.PSI(self.r(self.rnm, j, i)) - self.Alpha_n * self.PSI(self.r(self.rnp, j, i))) * (self.FiN(j, i) - self.Fi[j, i])
                st = 0.5 * self.F[self.S] * (self.Alpha_s * self.PSI(self.r(self.rsp, j, i)) - (1 - self.Alpha_s) * self.PSI(self.r(self.rsm, j, i))) * (self.Fi[j, i] - self.FiS(j, i))
                self.SuDC[j, i] = et + wt + nt + st



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
            print "SuDC coefficients:"
            print self.SuDC

    #end calculateTVDcoefficients


    def calculateTDMACoefficients(self, i):
        '''
        book pag 220
        Apply TDMA S to N sweeping W to E
        The discretisation equation is given by
       
        In the book they have it reversed "j" is for lines and "i" for columns 
        '''
        #calculate on each vertical  from S -> N 
        for j in range(0, self.Ny):

            #Compute the TDMA coefficients
            self.alp[j] = self.aN[j, i].copy()
            self.bet[j] = self.aS[j, i].copy()
            self.D[j] = self.aP[j, i].copy()

            #the free term
            Su = self.SuW[j, i] + self.SuE[j, i] + self.SuS[j, i] + self.SuN[j, i] + self.SuDC[j, i]
            #Avoid problems at boundaries by calling a function which considers the boundary limitation on index

            #boundary conditions are set through the term C[j]   
            self.C[j] = self.aW[j, i] * self.FiW(j, i) + self.aE[j, i] * self.FiE(j, i) + Su
        #end for j


    #end calculateTDMACoefficients

    def callTVD(self):
        #solve with TDMA for this i vertical line S-N
        n = self.D.size
        iter = 0



        #Because we don't know the values of Fi in the middle the first calculation will be far off
        #Therefore we set an iterative cycle to calculate the values of Fi  
        while  self.maxiter > iter :

            self.calculateTVDcoefficients()
            #Swipe  from W to E
            for i in range(0, self.Nx - 1):

                #calculate the TDMA coefficients for column i
                self.calculateTDMACoefficients(i)

                if self.debug == True:
                    print "beta:", self.bet
                    print "D", self.D
                    print "alp", self.alp
                    print "C", self.C


                x = thomas(n, -self.bet[1:], self.D, -self.alp[:-1], self.C)
                self.Fi[:, i] = x.copy()
            #end i

            #test accuracy and exit condition 
            flat = self.Fi[1] - self.FiOld[1]
            dx = math.sqrt(numpy.dot(flat, flat))

            print "iteration # %d, dx=%f" % (iter, dx)
            #print "Fi:", self.Fi

            #Exit if we are satisfied wit the accuracy
            if dx < self.EPSILON :
                return

            iter += 1

            #copy current values to the old values matrix
            self.FiOld = self.Fi.copy()

            #if we did not converge yet print an error and exit
            if self.maxiter < iter:
                print "Max iterations exceeded => did not converge"
                return
        #end while 
    #end callTDMA


    def plotFi(self):

        #create the figure and axes object
        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
        ax = fig.add_subplot(111)
#
        # Creating the grid of coordinates x,y 
        x = numpy.arange(0., self.Lx, self.deltaX)              #create the x, and y divisions
        y = numpy.arange(0., self.Ly, self.deltaY)
        X, Y = numpy.meshgrid(x, y)                             #create the grid


        #plot a heatmap
        im = ax.pcolor(X, Y, self.Fi)                           # draw the heatmap
        fig.colorbar(im)                                        # add the legend on a colour bar

        #superimpose contours
        cs = ax.contour(X, Y, self.Fi)                          # draw the contour 
        plt.clabel(cs, inline = 1, fontsize = 10)               # draw the contour

        #draw only the contours
        fig2 = plt.figure(2, facecolor = 'w', edgecolor = 'k')
        ax2 = fig2.add_subplot(111)
        cs = ax2.contour(X, Y, self.Fi)                          # draw the contour
        plt.clabel(cs, inline = 1, fontsize = 10)                # create the label on the contour


        plt.show()
    #end plotFi

    def solve(self):

        #We know E,W, N, S boundary conditions, so lets speed up things
        self.setBoundaryConditions()


        self.callTVD()

        #print the solution vector
        print self.Fi


    #end solve

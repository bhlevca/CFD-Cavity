'''
Created on Feb 6, 2012

@author: bogdan hlevca 995151213
'''

#import numerical Python for Matrix operations
import numpy
import SIMPLE
import linalg

#import the TDMA module
from thomas import *

#import the graphi library
import matplotlib.pyplot as plt

# Create a mesh class that holds a vector of nodes


#Class definition
class FinVol_2D_Conv_Diff(object):
    '''
    Solve the equation:
    d/dx(rho*u*Fi) + d/dy(ro*v*Fi) = d/dx(Gamma* dFi/dx)+d/dy(Gamma* dFi/dy) 
    
    where 
        d is the partial difference d_rond
        F - stands for the property Fi
    '''

    @classmethod
    def plot3curves(cls, FiCD, FiUD, FiQU):
        '''
        Plots comparative values on the diagonal of the square
        '''
        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
        ax = fig.add_subplot(111)

        #create the figure and axes object

        # Creating the grid of coordinates x,y 
        x = range(0, FiCD.size)                #create the x, and y divisions

        #draw only the contours
        cs = ax.plot(x, FiCD, x, FiUD, x, FiQU)                # plot diagonal values
        ax.legend(('Central Differences', 'Upwind Differences', 'QUICK'))
        plt.show()
    #end plot3curves

    #constructor object initializes the class arrays and calculates preliminary coefficients
    def __init__(self, var, FiBound, p, simple):
        '''
        Constructor
        '''

        self.debug = simple.debug
        self.simple = simple
        self.var = var
        self.W = 0
        self.E = 1
        self.N = 2
        self.S = 3

        self.urf = 0.7  #underrelaxation factor

        self.Gamma = numpy.zeros(4)
        self.A = numpy.zeros(4)

        self.Df = numpy.zeros(4)

        self.deltaX = None
        self.deltaY = None

        #iterative error accepted
        self.EPSILON = 1e-7
        self.maxiter = 10000

        self.Fi0E = FiBound[self.E]
        self.Fi0N = FiBound[self.N]
        self.Fi0W = FiBound[self.W]
        self.Fi0S = FiBound[self.S]

        self.Lx = simple.Lx
        self.Ly = simple.Ly
        self.Nx = simple.Nx
        self.Ny = simple.Ny
        self.p = p
        self.rho = simple.rho
        self.scheme = simple.scheme


        self.initNodes()               #initialize node coefficients
        self.calculateGamma(simple.Gx, simple.Gy)    #calculate preliminary information, Gamma, deltaX and Y and the control volume area 
        self.calculateDelta()
        self.calculateA()
        # convection terms are provided from outside by simple
        # self.calculate_F()
        self.calculate_Df()

        if self.simple.lc.CPP == True:
            self.cls = linalg.CppBlas()

    #end __init__

    def initNodes(self):
        '''
        initialize node coefficients and TDMA coefficients
        '''
        #init the nodes


        #must have their own coefficients as they differ for various variables
        self.aW = numpy.zeros((self.Ny, self.Nx))
        self.aE = numpy.zeros((self.Ny, self.Nx))
        self.aN = numpy.zeros((self.Ny, self.Nx))
        self.aS = numpy.zeros((self.Ny, self.Nx))
        self.aP = numpy.zeros((self.Ny, self.Nx))

        self.aWW = numpy.zeros((self.Ny, self.Nx))
        self.aEE = numpy.zeros((self.Ny, self.Nx))
        self.aNN = numpy.zeros((self.Ny, self.Nx))
        self.aSS = numpy.zeros((self.Ny, self.Nx))

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


        #init with zero the solution matrix     
        self.Fi = numpy.zeros((self.Ny, self.Nx))
        self.FiOld = numpy.zeros((self.Ny, self.Nx))

        #initialize with values from previous run if any
        if self.var == "u":
            self.Fi = self.simple.ustar.copy()

        elif self.var == "v":
            self.Fi = self.simple.vstar.copy()
        else:
            print "Wrong variable passed"
            raise NameError("Wrong variable passed")

        #set the old values as current
        self.FiOld = self.Fi.copy()

        #convection terms
        self.F = self.simple.F
        self.FOld = self.simple.FOld


        #the TDMA coeficients. (We calculate vertical lines (columns = > allocate the number of horizontal lines)
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
        if self.Nx == 1:
            self.A[self.W] = self.A[self.E] = 1
            self.A[self.N] = self.A[self.S] = 1
        else:
            self.A[self.W] = self.A[self.E] = self.deltaY
            self.A[self.N] = self.A[self.S] = self.deltaX

    def calculateDelta(self):
        '''
        calculate preliminary information, deltaX and Y 
        '''
        self.deltaX = self.Lx / self.Nx
        self.deltaY = self.Ly / self.Ny

    def calculate_F(self):
        '''unsteady flow F is NOT constant
       
        '''
        self.simple.calculate_F()

    #end Calculate_F



    def calculate_Df(self):
        '''
        calculate the diffusion coeff D
        '''
        self.Df[self.W] = self.Gamma[self.W] * self.A[self.W] / self.deltaX
        self.Df[self.E] = self.Gamma[self.E] * self.A[self.E] / self.deltaX
        self.Df[self.S] = self.Gamma[self.S] * self.A[self.S] / self.deltaY
        self.Df[self.N] = self.Gamma[self.N] * self.A[self.N] / self.deltaY

        print "***********************"
        print "Df on boundary set to 0"
        print "***********************"

        if self.scheme == "QUICK":
            self.D0W = self.Df[self.W]
            self.D0E = self.Df[self.E]
            self.D0S = self.Df[self.S]
            self.D0N = self.Df[self.N]
        else:
            self.D0W = 2 * self.Df[self.W]
            self.D0E = 2 * self.Df[self.E]
            self.D0S = 2 * self.Df[self.S]
            self.D0N = 2 * self.Df[self.N]

    #end calculate_D

    def calculate_alpha(self, j, i):
        if self.F[self.W, j, i] > 0:
            self.Alpha_w = 1
        else :
            self.Alpha_w = 0

        if self.F[self.E, j, i] > 0:
            self.Alpha_e = 1
        else :
            self.Alpha_e = 0

        if self.F[self.S, j, i] > 0:
            self.Alpha_s = 1
        else :
            self.Alpha_s = 0

        if self.F[self.N, j, i] > 0:
            self.Alpha_n = 1
        else :
            self.Alpha_n = 0
    #end calculate_alpha



    def setBoundaryConditions(self):
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                #west
                if i == 0:
                    self.Fi[j, i] = self.Fi0W

                #east
                if i == self.Nx - 1:
                    self.Fi[j, i] = self.Fi0E

                #north
                if j == self.Ny - 1:
                    self.Fi[j, i] = self.Fi0N

                #south
                if j == 0:
                    self.Fi[j, i] = self.Fi0S


    #end set BoundaryConditions


    # set of functions to account for the cardinal values of the main feature Fi
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
            FiWW = 0#2 * self.Fi0W - self.Fi[j, i]
        elif i == 1:
            FiWW = self.Fi0W
        else:
            FiWW = self.Fi[j, i - 2]
        return FiWW

    def FiEE(self, j, i):
        if i == self.Nx - 1  :
            FiEE = 0 #2 * self.Fi0E - self.Fi[j, i]
        elif  i == self.Nx - 2:
            FiEE = self.Fi0E
        else:
            FiEE = self.Fi[j, i + 2]
        return FiEE

    def FiSS(self, j, i):
        if j == 0:
            FiSS = 0 #2 * self.Fi0S - self.Fi[j, i]
        elif j == 1:
            FiSS = self.Fi0S
        else:
            FiSS = self.Fi[j - 2, i]
        return FiSS

    def FiNN(self, j, i):
        if j == self.Ny - 1  :
            FiNN = 0 # 2 * self.Fi0N - self.Fi[j, i]
        elif j == self.Ny - 2:
            FiNN = self.Fi0N
        else:
            FiNN = self.Fi[j + 2, i]
        return FiNN




    def calculate_C(self, j, i, aW, FiW, aE, FiE, Su, p , var):
        if (var == "u"):
            return aW[j, i] * FiW(j, i) + aE[j, i] * FiE(j, i) + (self.simple.pw(p, j, i) - self.simple.pe(p, j, i)) * self.A[self.E] + Su
        elif (var == "v"):
            return aW[j, i] * FiW(j, i) + aE[j, i] * FiE(j, i) + (self.simple.ps(p, j, i) - self.simple.pn(p, j, i)) * self.A[self.N] + Su


    def calculateTDMACoefficients(self, i):
        '''
        book pag 220
        Apply TDMA S to N sweeping W to E
        The discretization equation is given by
       
        In the book they have it reversed "j" is for lines and "i" for columns 
        '''
        #calculate on each vertical  from S -> N 

        for j in range(0, self.Ny):

            #Compute the TDMA coefficients
            self.alp[j] = self.aN[j, i].copy()
            self.bet[j] = self.aS[j, i].copy()
            self.D[j] = self.aP[j, i].copy()

            #the free term
            #Avoid problems at boundaries by calling a function which considers the boundary limitation on index

            #boundary conditions are set through the term C[j] 

            Su = self.SuW[j, i] + self.SuE[j, i] + self.SuS[j, i] + self.SuN[j, i]\
                + self.aWW[j, i] * self.FiWW(j, i) + self.aEE[j, i] * self.FiEE(j, i) + self.aNN[j, i] * self.FiNN(j, i) + self.aSS[j, i] * self.FiSS(j, i)

            #self.C[j] = self.aW[j, i] * self.FiW(j, i) + self.aE[j, i] * self.FiE(j, i) + (self.pE(j,i) - self.pE(j,i)) + Su
            self.C[j] = self.calculate_C(j, i, self.aW, self.FiW, self.aE, self.FiE, Su, self.p, self.var)
        #end for j


    #end calculateTDMACoefficients



    def calculateQUICKSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions for QUICK scheme
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #self.calculate_alpha(j, i)
                #west
                if i == 0:
                    #cof = (8.0 / 3 * self.Df[self.W] + 2.0 / 8 * self.F[self.E, j, i] + self.F[self.W, j, i])
                    cof = (8.0 / 3 * self.D0W + 2.0 / 8 * self.F[self.E, j, i] + self.F[self.W, j, i])
                    self.SuW[j, i] = cof * self.Fi0W
                    self.SpW[j, i] = -cof
                elif i == 1 :
                    cof = 1.0 / 4 * self.F[self.W, j, i]
                    self.SuW[j, i] = -cof * self.Fi0W
                    self.SpW[j, i] = cof
                else:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0

                #east
                if i == self.Nx - 1:
                    cof = (8.0 / 3 * self.D0E - self.F[self.E, j, i])
                    self.SuE[j, i] = cof * self.Fi0E
                    self.SpE[j, i] = -cof
                else:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0

                #south
                if j == 0:
                    cof = (8.0 / 3 * self.D0S + 2.0 / 8 * self.F[self.S, j, i] + self.F[self.S, j, i])
                    self.SuS[j, i] = cof * self.Fi0S
                    self.SpS[j, i] = -cof
                elif j == 1:
                    cof = 1.0 / 4 * self.F[self.S, j, i]
                    self.SuS[j, i] = -cof * self.Fi0S
                    self.SpS[j, i] = cof
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0

                #north
                if j == self.Ny - 1:
                    cof = (8.0 / 3 * self.D0N - self.F[self.N, j, i])
                    self.SuN[j, i] = cof * self.Fi0N
                    self.SpN[j, i] = -cof
                else:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0

        #end j
        if self.debug == True:
            if self.Nx == 1:
                print "Su:", self.SuS + self.SuN
                print "Sp:", self.SpS + self.SpN

            else:
                print "Su:", self.SuE + self.SuW + self.SuS + self.SuN
                print "Sp:", self.SpE + self.SpW + self.SpS + self.SpN

    #end calculateStandardSources



    def calculateCDSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions for central differences scheme
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #west
                if i == 0:
                    self.SuW[j, i] = (self.D0W + self.F[self.W, j, i]) * self.Fi0W
                    self.SpW[j, i] = -(self.D0W + self.F[self.W, j, i])
                else:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0

                #east
                if i == self.Nx - 1:
                    self.SuE[j, i] = (self.D0E - self.F[self.E, j, i]) * self.Fi0E
                    self.SpE[j, i] = -(self.D0E - self.F[self.E, j, i])
                else:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0

                #north 
                if j == self.Ny - 1:
                    self.SuN[j, i] = (self.D0N - self.F[self.N, j, i]) * self.Fi0N
                    self.SpN[j, i] = -(self.D0N - self.F[self.N, j, i])
                else:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0
                #south
                if j == 0:
                    self.SuS[j, i] = (self.D0S + self.F[self.S, j, i]) * self.Fi0S
                    self.SpS[j, i] = -(self.D0E + self.F[self.S, j, i])
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
        #end j
        if self.debug == True:
            print "Su:", self.SuE + self.SuW + self.SuS + self.SuN
            print "Sp:", self.SpE + self.SpW + self.SpS + self.SpN

    #end calculateCDSources

    def calculateUDSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions for upwind differences scheme
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #west
                if i == 0:
                    #self.SuW[j, i] = (2 * self.Df[self.W] + self.F[self.W, j, i]) * self.Fi0W
                    #self.SpW[j, i] = -(2 * self.Df[self.W] + self.F[self.W, j, i])

                    cof = (self.D0W + self.F[self.W, j, i])
                    self.SuW[j, i] = cof * self.Fi0W
                    self.SpW[j, i] = -cof
                else:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0

                #east
                if i == self.Nx - 1:
                    cof = self.D0E
                    self.SuE[j, i] = cof * self.Fi0E
                    self.SpE[j, i] = -cof
                    #self.SuE[j, i] = 2 * self.Df[self.E] * self.Fi0E
                    #self.SpE[j, i] = -2 * self.Df[self.E]
                else:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0

                #north 
                if j == self.Ny - 1:
                    #self.SuN[j, i] = 2 * self.Df[self.N] * self.Fi0N
                    #self.SpN[j, i] = -2 * self.Df[self.N]
                    cof = self.D0N
                    self.SuN[j, i] = cof * self.Fi0N
                    self.SpN[j, i] = -cof
                else:
                    self.SuN[j, i] = 0
                    self.SpN[j, i] = 0
                #south
                if j == 0:
                    #self.SuS[j, i] = (2 * self.Df[self.S] + self.F[self.S, j, i]) * self.Fi0S
                    #self.SpS[j, i] = -(2 * self.Df[self.S] + self.F[self.S, j, i])
                    cof = self.D0S + self.F[self.S, j, i]
                    self.SuS[j, i] = cof * self.Fi0S
                    self.SpS[j, i] = -cof
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0
        #end j
        if self.debug == True:
            print "Su:", self.SuE + self.SuW + self.SuS + self.SuN
            print "Sp:", self.SpE + self.SpW + self.SpS + self.SpN

    #end calculateuDSources


    def calculate_aP(self, scheme, Sp, fE, fW, fN, fS, aE, aW, aN, aS, aEE = 0, aWW = 0, aNN = 0, aSS = 0):
        if scheme == 'QUICK':
            return aW + aE + aS + aN + (fE - fW) + (fN - fS) - Sp + aSS + aNN + aWW + aEE
        elif scheme == 'UD':
            return aW + aE + aS + aN + (fE - fW) + (fN - fS) - Sp
        elif scheme == 'CD':
            return aW + aE + aS + aN + (fE - fW) + (fN - fS) - Sp
    #end calculate_aP


    def calculateQUICKCoefficients(self):
        '''
            2D "a" coeficients for QUICK are implementation at page 163 in the Versteeg book 
            
        '''
        self.calculateQUICKSources()

        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                self.calculate_alpha(j, i)
                #aW
                if i == 0:
                    self.aW[j, i] = 0
                elif i == 1:
                    self.aW[j, i] = self.Df[self.W] + 7.0 / 8 * self.F[self.W, j, i] + 1.0 / 8 * self.F[self.E, j, i]
                elif i == self.Nx - 1:
                    self.aW[j, i] = self.Df[self.W] + 1.0 / 3 * self.D0E + 6.0 / 8 * self.F[self.W, j, i]
                else:
                    self.aW[j, i] = self.Df[self.W] + 6.0 / 8 * self.Alpha_w * self.F[self.W, j, i] + 1.0 / 8 * self.Alpha_e * self.F[self.E, j, i]\
                     + 3.0 / 8 * (1 - self.Alpha_w) * self.F[self.W, j, i]

                #aE
                if i == 0 :
                    self.aE[j, i] = self.Df[self.E] + 1.0 / 3 * self.D0W - 3.0 / 8 * self.F[self.E, j, i]
                elif i == 1:
                    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E, j, i]
                elif i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E, j, i] - 6.0 / 8 * (1 - self.Alpha_e) * self.F[self.E, j, i]\
                     - 1.0 / 8 * (1 - self.Alpha_w) * self.F[self.W, j, i]

                #aWW
                if i == 0:
                    self.aWW[j, i] = 0
                elif i == 1:
                    self.aWW[j, i] = 0
                else:
                    self.aWW[j, i] = -1.0 / 8 * self.Alpha_w * self.F[self.W, j, i]

                #aEE
                if i == self.Nx - 1 or i == self.Nx - 2:
                    self.aEE[j, i] = 0
                else:
                    self.aEE[j, i] = 1.0 / 8 * (1 - self.Alpha_e) * self.F[self.E, j, i]

                #aS
                if j == 0:
                    self.aS[j, i] = 0
                elif j == 1:
                    self.aS[j, i] = self.Df[self.S] + 7.0 / 8 * self.F[self.S, j, i] + 1.0 / 8 * self.F[self.N, j, i]
                elif j == self.Ny - 1:
                    self.aS[j, i] = self.Df[self.S] + 1.0 / 3 * self.D0N + 6.0 / 8 * self.F[self.S, j, i]
                else:
                    self.aS[j, i] = self.Df[self.S] + 6.0 / 8 * self.Alpha_s * self.F[self.S, j, i] + 1.0 / 8 * self.Alpha_n * self.F[self.N, j, i] + 3.0 / 8 * (1 - self.Alpha_s) * self.F[self.S, j, i]

                #aNself.eps
                if j == 0 :
                    self.aN[j, i] = self.Df[self.N] + 1.0 / 3 * self.D0S - 3.0 / 8 * self.F[self.N, j, i]
                elif j == 1:
                    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N, j, i]
                elif j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N, j, i] - 6.0 / 8 * (1 - self.Alpha_n) * self.F[self.N, j, i] - 1.0 / 8 * (1 - self.Alpha_s) * self.F[self.S, j, i]

                #aSS
                if j == 0:
                    self.aSS[j, i] = 0
                elif j == 1:
                    self.aSS[j, i] = 0
                else:
                    self.aSS[j, i] = -1.0 / 8 * self.Alpha_s * self.F[self.S, j, i]

                #aNN
                if j == self.Ny - 1 or j == self.Ny - 2:
                    self.aNN[j, i] = 0
                else:
                    self.aNN[j, i] = 1.0 / 8 * (1 - self.Alpha_n) * self.F[self.N, j, i]


                Sp = self.SpE[j, i] + self.SpW[j, i] + self.SpS[j, i] + self.SpN[j, i]
                #self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + self.aSS[j, i] + self.aNN[j, i] + self.aWW[j, i] + self.aEE[j, i]\
                #     + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp

                self.aP[j, i] = self.calculate_aP('QUICK', Sp, \
                               self.F[self.E, j, i], self.F[self.W, j, i], self.F[self.N, j, i], self.F[self.S, j, i], \
                               self.aE[j, i], self.aW[j, i], self.aN[j, i], self.aS[j, i], \
                               self.aEE[j, i], self.aWW[j, i], self.aNN[j, i], self.aSS[j, i])

                #TVD deferred correction source term
            #end i
        #end j

    #end calculateStandardOUICKcoefficients


    def calculateUDCoefficients(self):
        '''
        calculare matrix coefficients for Upwind scheme
        '''
        self.calculateUDSources()

        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                #west
                if i == 0:
                    self.aW[j, i] = 0
                else:
                    self.aW[j, i] = self.Df[self.W] + max(self.F[self.W, j, i], 0)

                #east
                if i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    self.aE[j, i] = self.Df[self.E] + max(0, -self.F[self.E, j, i])

                #north 
                if j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    self.aN[j, i] = self.Df[self.N] + max(0, -self.F[self.N, j, i])
                #south
                if j == 0:
                    self.aS[j, i] = 0
                else:
                    self.aS[j, i] = self.Df[self.S] + +max(self.F[self.S, j, i], 0)

                #ap coefficient
                Sp = self.SpW[j, i] + self.SpE[j, i] + self.SpS[j, i] + self.SpN[j, i]
                #self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp
                self.aP[j, i] = self.calculate_aP('UD', Sp, \
                               self.F[self.E, j, i], self.F[self.W, j, i], self.F[self.N, j, i], self.F[self.S, j, i], \
                               self.aE[j, i], self.aW[j, i], self.aN[j, i], self.aS[j, i])
        #END for i
    #end calculateUDCoefficients 

    def calculateCDCoefficients(self):
        '''
        calculate matrix coeff for upwinding scheme
        '''
        self.calculateCDSources()

        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                #west
                if i == 0:
                    self.aW[j, i] = 0
                else:
                    self.aW[j, i] = self.Df[self.W] + self.F[self.W, j, i] / 2

                #east
                if i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    self.aE[j, i] = self.Df[self.E] - self.F[self.E, j, i] / 2

                #north 
                if j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    self.aN[j, i] = self.Df[self.N] - self.F[self.N, j, i] / 2
                #south
                if j == 0:
                    self.aS[j, i] = 0
                else:
                    self.aS[j, i] = self.Df[self.S] + self.F[self.S, j, i] / 2

                #ap coefficient
                Sp = self.SpW[j, i] + self.SpE[j, i] + self.SpS[j, i] + self.SpN[j, i]
                #self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp
                self.aP[j, i] = self.calculate_aP('CD', Sp, \
                               self.F[self.E, j, i], self.F[self.W, j, i], self.F[self.N, j, i], self.F[self.S, j, i], \
                               self.aE[j, i], self.aW[j, i], self.aN[j, i], self.aS[j, i])
        #END for i
    #end calculateCDCoefficients

    def callSolver(self):
        #solve with TDMA for this i vertical line S-N
        n = self.D.size
        Iter = 0

        #does not depend on Fi so we take it our of the iteration process
        if self.scheme == "QUICK":
            self.calculateQUICKCoefficients()
        elif self.scheme == "CD":
            self.calculateCDCoefficients()
        elif self.scheme == "UD":
            self.calculateUDCoefficients()
        else:
            print "Unknown Scheme!!!"
            return

        #Because we don't know the values of Fi in the middle the first calculation will be far off
        #Therefore we set an iterative cycle to calculate the values of Fi  

        if self.simple.lc.CPP == True:
            x = numpy.zeros(n)

        while  self.maxiter > Iter :

            #Swipe  from W to E
            for i in range(0, self.Nx):

                #calculate the TDMA coefficients for column i
                self.calculateTDMACoefficients(i)

                if self.debug == True:
                    print "beta:", self.bet
                    print "D", self.D
                    print "alp", self.alp
                    print "C", self.C

                if self.simple.lc.CPP == True:
                    self.cls.setTDMA(-self.bet[1:], self.D, -self.alp[:-1], self.C, n)
                    d = self.cls.solveTDMA(x, n)
                    self.Fi[:, i] = d["solution"].copy()
                else :
                    x = thomas(n, -self.bet[1:], self.D, -self.alp[:-1], self.C)
                    self.Fi[:, i] = x.copy()

            #end i

            #TODO Under relaxation
            self.Fi = self.urf * self.Fi.copy() + self.FiOld.copy() * (1 - self.urf)

            #test accuracy and exit condition
            if self.Nx == 1:
                flat = self.Fi[:, 0] - self.FiOld[:, 0]
            else:
                flat = self.Fi[:, 1] - self.FiOld[:, 1]
            dx = math.sqrt(numpy.dot(flat, flat))

            if Iter % 600 == 0:
                print "var: %s iter # %d, dx=%1.9f" % (self.var, Iter, dx)
            #print "Fi:", self.Fi

            #Exit if we are satisfied wit the accuracy
            if dx < self.EPSILON :
                return

            Iter += 1

            #copy current values to the old values matrix
            self.FiOld = self.Fi.copy()

            #if we did not converge yet print an error and exit
            if self.maxiter < Iter:
                print "Max iterations exceeded => did not converge"
                return
        #end while 
    #end callQUICK


    def plotFi(self):

        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
        ax = fig.add_subplot(111)

        if self.Nx == 1:
            y = numpy.arange(0., self.Ly, self.deltaY)
            ax.plot(y, self.Fi[:, 0])
        else:
            #create the figure and axes object

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
        #end if

        plt.show()
    #end plotFi


    def solve(self):

        #We know E,W, N, S boundary conditions, so lets speed up things
        self.setBoundaryConditions()

        #call the main algorithm
        self.callSolver()

        #print the solution vector
        print self.Fi
        return self.Fi

    #end solve



#end class





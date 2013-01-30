'''
Created on Feb 6, 2012

@author: bogdan
'''

#import numerical Python for Matrix operations
import numpy

#import the TDMA module
from thomas import *

#import the graphi library
import matplotlib.pyplot as plt


#Class definition
class FinVol_QUICK_Hayase_2D_Conv_Diff:
    '''
    Solve the equation:
    d/dx(ro*u*F) + d/dy(ro*v*F) = d/dx(Gamma* dF/dx)+d/dy(Gamma* dF/dy) 
    
    where 
        d is the partial difference d_rond
        F - stands for the property Fi
    '''

    #constructor object initializes the class arrays and calculates preliminary coefficients
    def __init__(self, Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, debug = False):
        '''
        Constructor
        '''

        self.debug = debug
        self.W = 0
        self.E = 1
        self.N = 2
        self.S = 3


        self.Gamma = numpy.zeros(4)
        self.A = numpy.zeros(4)
        self.F = numpy.zeros(4)
        self.Df = numpy.zeros(4)

        self.deltaX = None
        self.deltaY = None

        #iterative error accepted
        self.EPSILON = 1e-7
        self.maxiter = 2000

        self.Fi0E = FiE
        self.Fi0N = FiN
        self.Fi0W = FiW
        self.Fi0S = FiS

        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
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
        #init Su deferred sources
        self.SuDC = numpy.zeros((self.Ny, self.Nx))

        #init with zero the solution matrix     
        self.Fi = numpy.zeros((self.Ny, self.Nx))
        self.FiOld = numpy.zeros((self.Ny, self.Nx))

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
        calculate the diffusion coeff D
        '''
        self.Df[self.W] = self.Gamma[self.W] * self.A[self.W] / self.deltaX
        self.Df[self.E] = self.Gamma[self.E] * self.A[self.E] / self.deltaX
        self.Df[self.S] = self.Gamma[self.S] * self.A[self.S] / self.deltaY
        self.Df[self.N] = self.Gamma[self.N] * self.A[self.N] / self.deltaY
        self.D0W = self.Df[self.W]
        self.D0E = self.Df[self.E]
        self.D0S = self.Df[self.S]
        self.D0N = self.Df[self.N]
        
    #end calculate_Df

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


    
    def setBoundaryConditions(self):
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                #west
                if i == 0:
                    self.Fi[j, i] = self.Fi0W

                #east
                if i == self.Nx - 1:
                    self.Fi[j, i] = self.Fi0E

                #northself.Su
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
            FiWW = 2 * self.Fi0W - self.Fi[j, i]
        elif i == 1:
            FiWW = self.Fi0W
        else:
            FiWW = self.Fi[j, i - 2]
        return FiWW

    def FiEE(self, j, i):
        if i == self.Nx - 1  :
            FiEE = 2 * self.Fi0E - self.Fi[j, i]
        elif  i == self.Nx - 2:
            FiEE = self.Fi0E
        else:
            FiEE = self.Fi[j, i + 2]
        return FiEE

    def FiSS(self, j, i):
        if j == 0:
            FiSS = 2 * self.Fi0S - self.Fi[j, i]
        elif j == 1:
            FiSS = self.Fi0S
        else:
            FiSS = self.Fi[j - 2, i]
        return FiSS

    def FiNN(self, j, i):
        if j == self.Ny - 1  :
            FiNN = 2 * self.Fi0N - self.Fi[j, i]
        elif j == self.Ny - 2:
            FiNN = self.Fi0N
        else:
            FiNN = self.Fi[j + 2, i]
        return FiNN

    def calculateQUICKSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions for QUICK scheme
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                #west
                if i == 0:
                    cof = (8.0 / 3 * self.Df[self.W] + 2.0 / 8 * self.F[self.E] + self.F[self.W])
                    self.SuW[j, i] = cof * self.Fi0W
                    self.SpW[j, i] = -cof
                elif i == 1 :
                    cof = 1.0 / 4 * self.F[self.W]
                    self.SuW[j, i] = -cof * self.Fi0W
                    self.SpW[j, i] = cof
                else:
                    self.SuW[j, i] = 0
                    self.SpW[j, i] = 0

                #east
                if i == self.Nx - 1:
                    cof = (8.0 / 3 * self.Df[self.E] - self.F[self.E])
                    self.SuE[j, i] = cof * self.Fi0E
                    self.SpE[j, i] = -cof
                else:
                    self.SuE[j, i] = 0
                    self.SpE[j, i] = 0

                #south
                if j == 0:
                    cof = (8.0 / 3 * self.Df[self.S] + 2.0 / 8 * self.F[self.S] + self.F[self.S])
                    self.SuS[j, i] = cof * self.Fi0S
                    self.SpS[j, i] = -cof
                elif j == 1:
                    cof = 1.0 / 4 * self.F[self.S]
                    self.SuS[j, i] = -cof * self.Fi0S
                    self.SpS[j, i] = cof
                else:
                    self.SuS[j, i] = 0
                    self.SpS[j, i] = 0

                #north
                if j == self.Ny - 1:
                    cof = (8.0 / 3 * self.Df[self.N] - self.F[self.N])
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

    #--------------------------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------------------------
    def O_calculateQUICKCoefficients(self):
        '''
            2D "a" coeficients for QUICK are implementation at page 163 in the Versteeg book 
            
        '''
        self.calculateQUICKSources()


        for j in range(0, self.Ny):
            for i in range(0, self.Nx) :
                #aW
                if i == 0:
                    self.aW[j, i] = 0
                elif i == 1:
                    self.aW[j, i] = self.Df[self.W] + 7.0 / 8 * self.F[self.W] + 1.0 / 8 * self.F[self.E]
                elif i == self.Nx - 1:
                    self.aW[j, i] = self.Df[self.W] + 1.0 / 3 * self.D0E + 6.0 / 8 * self.F[self.W]
                else:
                    self.aW[j, i] = self.Df[self.W] + 6.0 / 8 * self.Alpha_w * self.F[self.W] + 1.0 / 8 * self.Alpha_e * self.F[self.E]\
                     + 3.0 / 8 * (1 - self.Alpha_w) * self.F[self.W]

                #aE
                if i == 0 :
                    self.aE[j, i] = self.Df[self.E] + 1.0 / 3 * self.D0W - 3.0 / 8 * self.F[self.E]
                elif i == 1:
                    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E]
                elif i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E] - 6.0 / 8 * (1 - self.Alpha_e) * self.F[self.E]\
                     - 1.0 / 8 * (1 - self.Alpha_w) * self.F[self.W]

                #aWW
                if i == 0:
                    self.aWW[j, i] = 0
                elif i == 1:
                    self.aWW[j, i] = 0
                else:
                    self.aWW[j, i] = -1.0 / 8 * self.Alpha_w * self.F[self.W]

                #aEE
                if i == self.Nx - 1 or i == self.Nx - 2:
                    self.aEE[j, i] = 0
                else:
                    self.aEE[j, i] = 1.0 / 8 * (1 - self.Alpha_e) * self.F[self.E]


                #aS
                if j == 0:
                    self.aS[j, i] = 0
                elif j == 1:
                    self.aS[j, i] = self.Df[self.S] + 7.0 / 8 * self.F[self.S] + 1.0 / 8 * self.F[self.N]
                elif j == self.Ny - 1:
                    self.aS[j, i] = self.Df[self.S] + 1.0 / 3 * self.D0N + 6.0 / 8 * self.F[self.S]
                else:
                    self.aS[j, i] = self.Df[self.S] + 6.0 / 8 * self.Alpha_s * self.F[self.S] + 1.0 / 8 * self.Alpha_n * self.F[self.N] + 3.0 / 8 * (1 - self.Alpha_s) * self.F[self.S]

                #aN
                if j == 0 :
                    self.aN[j, i] = self.Df[self.N] + 1.0 / 3 * self.D0S - 3.0 / 8 * self.F[self.N]
                elif j == 1:
                    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N]
                elif j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N] - 6.0 / 8 * (1 - self.Alpha_n) * self.F[self.N] - 1.0 / 8 * (1 - self.Alpha_s) * self.F[self.S]

                #aSS
                if j == 0:
                    self.aSS[j, i] = 0
                elif j == 1:
                    self.aSS[j, i] = 0
                else:
                    self.aSS[j, i] = -1.0 / 8 * self.Alpha_s * self.F[self.S]

                #aNN
                if j == self.Ny - 1 or j == self.Ny - 2:
                    self.aNN[j, i] = 0
                else:
                    self.aNN[j, i] = 1.0 / 8 * (1 - self.Alpha_n) * self.F[self.N]


                if self.Nx == 1:
                    Sp = self.SpS[j, i] + self.SpN[j, i]
                    self.aP[j, i] = self.aS[j, i] + self.aN[j, i] + self.aSS[j, i] + self.aNN[j, i] + (self.F[self.N] - self.F[self.S]) - Sp
                else :
                    Sp = self.SpE[j, i] + self.SpW[j, i] + self.SpS[j, i] + self.SpN[j, i]
                    self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + self.aSS[j, i] + self.aNN[j, i] + self.aWW[j, i] + self.aEE[j, i]\
                     + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp

                #TVD deferred correction source term
            #end i
        #end j
        if self.debug == True:
            print "aW coefficients:"
            print self.aW
            print "aE coefficients:"
            print self.aE
            print "aS coefficients:"
            print self.aS
            print "aN coefficients:"
            print self.aN

            print "aWW coefficients:"
            print self.aWW
            print "aEE coefficients:"
            print self.aEE
            print "aSS coefficients:"
            print self.aSS
            print "aNN coefficients:"
            print self.aNN


            print "aP coefficients:"
            print self.aP
    #end calculateStandardOUICKcoefficients            


    def calculateQUICKCoefficients(self):
        '''
            2D "a" coeficients for TVD are implementation at page 175 in the Versteeg book 
            
        '''
        self.calculateQUICKSources()
        
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):

                #aW
                if i == 0:
                    self.aW[j, i] = 0
                elif i == 1:
                    self.aW[j, i] = self.Df[self.W] + 7.0 / 8 * self.F[self.W] + 1.0 / 8 * self.F[self.E]
                elif i == self.Nx - 1:
                    self.aW[j, i] = self.Df[self.W] + 1.0 / 3 * self.D0E + 6.0 / 8 * self.F[self.W]
                #elif i == self.Nx - 2:    
                #    self.aW[j, i] = self.Df[self.W] + 6.0 / 8 * self.Alpha_w * self.F[self.W] + 1.0 / 8 * self.Alpha_e * self.F[self.E]\
                #     + 3.0 / 8 * (1 - self.Alpha_w) * self.F[self.W]
                else:
                    self.aW[j, i] = self.Df[self.W] + self.Alpha_w * self.F[self.W]
                #aE
                if i == 0 :
                    self.aE[j, i] = self.Df[self.E] + 1.0 / 3 * self.D0W - 3.0 / 8 * self.F[self.E]
                elif i == 1:
                    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E]
                elif i == self.Nx - 1:
                    self.aE[j, i] = 0
                #elif i == self.Nx - 2:
                #    self.aE[j, i] = self.Df[self.E] - 3.0 / 8 * self.Alpha_e * self.F[self.E] - 6.0 / 8 * (1 - self.Alpha_e) * self.F[self.E]\
                #     - 1.0 / 8 * (1 - self.Alpha_w) * self.F[self.W]
                else:
                    self.aE[j, i] = self.Df[self.E] - (1 - self.Alpha_e) * self.F[self.E]

                #aWW
                if i == 0:
                    self.aWW[j, i] = 0
                elif i == 1:
                    self.aWW[j, i] = 0
                else:
                    self.aWW[j, i] = -1.0 / 8 * self.Alpha_w * self.F[self.W]

                #aEE
                if i == self.Nx - 1 or i == self.Nx - 2:
                    self.aEE[j, i] = 0
                else:
                    self.aEE[j, i] = 1.0 / 8 * (1 - self.Alpha_e) * self.F[self.E]
                    
                    
                #aS
                if j == 0:
                    self.aS[j, i] = 0
                elif j == 1:
                    self.aS[j, i] = self.Df[self.S] + 7.0 / 8 * self.F[self.S] + 1.0 / 8 * self.F[self.N]
                elif j == self.Ny - 1:
                    self.aS[j, i] = self.Df[self.S] + 1.0 / 3 * self.D0N + 6.0 / 8 * self.F[self.S]
                #elif j == self.Ny - 2:
                #    self.aS[j, i] = self.Df[self.S] + 6.0 / 8 * self.Alpha_s * self.F[self.S] + 1.0 / 8 * self.Alpha_n * self.F[self.N] + 3.0 / 8 * (1 - self.Alpha_s) * self.F[self.S]
                else:
                    self.aS[j, i] = self.Df[self.S] + self.Alpha_s * self.F[self.S]
                
                
                #aN
                if j == 0 :
                    self.aN[j, i] = self.Df[self.N] + 1.0 / 3 * self.D0S - 3.0 / 8 * self.F[self.N]
                elif j == 1:
                    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N]
                elif j == self.Ny - 1:
                    self.aN[j, i] = 0
                #elif j == self.Ny - 2:    
                #    self.aN[j, i] = self.Df[self.N] - 3.0 / 8 * self.Alpha_n * self.F[self.N] - 6.0 / 8 * (1 - self.Alpha_n) * self.F[self.N] - 1.0 / 8 * (1 - self.Alpha_s) * self.F[self.S]
                else:
                    self.aN[j, i] = self.Df[self.N] - (1 - self.Alpha_n) * self.F[self.N]

                #aSS
                if j == 0:
                    self.aSS[j, i] = 0
                elif j == 1:
                    self.aSS[j, i] = 0
                else:
                    self.aSS[j, i] = -1.0 / 8 * self.Alpha_s * self.F[self.S]

                #aNN
                if j == self.Ny - 1 or j == self.Ny - 2:
                    self.aNN[j, i] = 0
                else:
                    self.aNN[j, i] = 1.0 / 8 * (1 - self.Alpha_n) * self.F[self.N]
                    
                Sp = self.SpE[j, i] + self.SpW[j, i] + self.SpS[j, i] + self.SpN[j, i]
                
                if i==0 or j==0 or i==1 or j==1 or i== self.Nx-1 or j==self.Ny-1 or i== self.Nx-2:# or j==self.Ny-2:
                    self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i] + self.aSS[j, i] + self.aNN[j, i] + self.aWW[j, i] + self.aEE[j, i]\
                     + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S]) - Sp
                else:     
                    self.aP[j, i] = self.aW[j, i] + self.aE[j, i] + self.aS[j, i] + self.aN[j, i]\
                              + (self.F[self.E] - self.F[self.W]) + (self.F[self.N] - self.F[self.S])

    #end calculateOUICKcoefficients

    

    def calculateDeferredSources(self):
        '''
        Determine sources on all cardinal points based on boundary conditions
        
        For the Assignment 3 all the sources  = 0 since there is no flux on the borders
        '''

        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                _SuEp = 1.0 / 8.0 * (self.FiW(j, i) + 2.0 * self.Fi[j, i] - 3.0 * self.FiE(j, i)) * self.Alpha_e * self.F[self.E] 
                _SuWp = 1.0 / 8.0 * (3.0 * self.Fi[j, i] - 2.0 * self.FiW(j, i) - self.FiWW(j, i)) * self.Alpha_w * self.F[self.W]
                _SuWm = 1.0 / 8.0 * (3.0 * self.FiW(j, i) - 2.0 * self.Fi[j, i] - self.FiE(j, i)) * (1.0 - self.Alpha_w) * self.F[self.W]
                _SuEm = 1.0 / 8.0 * (2.0 * self.FiE(j, i) + self.FiEE(j, i) - 3.0 * self.Fi[j, i]) * (1.0 - self.Alpha_e) * self.F[self.E]

                _SuSp = 1.0 / 8.0 * (3.0 * self.Fi[j, i] - 2.0 * self.FiS(j, i) - self.FiSS(j, i)) * self.Alpha_s * self.F[self.S]
                _SuNp = 1.0 / 8.0 * (self.FiS(j, i) + 2.0 * self.Fi[j, i] - 3.0 * self.FiN(j, i)) * self.Alpha_n * self.F[self.N]
                _SuSm = 1.0 / 8.0 * (3.0 * self.FiS(j, i) - 2.0 * self.Fi[j, i] - self.FiN(j, i)) * (1.0 - self.Alpha_s) * self.F[self.S]
                _SuNm = 1.0 / 8.0 * (2.0 * self.FiN(j, i) + self.FiNN(j, i) - 3.0 * self.Fi[j, i]) * (1.0 - self.Alpha_n) * self.F[self.N]

                self.SuDC[j, i] = _SuWp + _SuEp + _SuEm + _SuWm + _SuSp + _SuNp + _SuSm + _SuNm
            #end j
        #end i

    #end calculateSources 

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
            #Avoid problems at boundaries by calling a function which considers the boundary limitation on index

            #boundary conditions are set through the term C[j] 
            
            if self.Nx == 1:
                self.C[j] = self.SuDC[j, i]
            else :
                if i == 0 or i==1 or j==0 or j == 1 or i == self.Nx-1 or j == self.Ny-1:# or i == self.Nx-2 or j == self.Ny-2:
                    Su = self.SuW[j, i] + self.SuE[j, i] + self.SuS[j, i] + self.SuN[j, i]\
                        + self.aWW[j, i] * self.FiWW(j, i) + self.aEE[j, i] * self.FiEE(j, i)\
                        + self.aNN[j, i] * self.FiNN(j, i) + self.aSS[j, i] * self.FiSS(j, i)
                else:
                    Su = self.SuDC[j, i]    
                
                self.C[j] = self.aW[j, i] * self.FiW(j, i) + self.aE[j, i] * self.FiE(j, i) + Su
                
        #end for j


    #end calculateTDMACoefficients

    def callTDMA(self):
        #solve with TDMA for this i vertical line S-N
        n = self.D.size
        Iter = 0

        #does not depend on Fi so we take it our of the iteration process
        self.calculateQUICKCoefficients()
        if self.debug == True:
            print "aW coefficients:"
            print self.aW
            print "aE coefficients:"
            print self.aE
            print "aS coefficients:"
            print self.aS
            print "aN coefficients:"
            print self.aN

            print "aWW coefficients:"
            print self.aWW
            print "aEE coefficients:"
            print self.aEE
            print "aSS coefficients:"
            print self.aSS
            print "aNN coefficients:"
            print self.aNN


            print "aP coefficients:"
            print self.aP

        
        
        #Because we don't know the values of Fi in the middle the first calculation will be far off
        #Therefore we set an iterative cycle to calculate the values of Fi  
        while  self.maxiter > Iter :

            
            self.calculateDeferredSources()

            if self.debug == True:
                print "SuDC coefficients:"
                print self.SuDC
            
            #Swipe  from W to E
            for i in range(0, self.Nx):

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
            if self.Nx == 1:
                flat = self.Fi[:, 0] - self.FiOld[:, 0]
            else:
                flat = self.Fi[:, 1] - self.FiOld[:, 1]
            dx = math.sqrt(numpy.dot(flat, flat))

            print "iteration # %d, dx=%1.9f" % (Iter, dx)
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
    #end callTDMA


    def plotFi(self):

        fig = plt.figure(1, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
        ax = fig.add_subplot(111)

        if self.Nx == 1:
            y = numpy.arange(0., self.Ly, self.deltaY)
            ax.plot(y, self.Fi[:, 0])
        else:
            #create the figure and axes object

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
        #end if

        plt.show()
    #end plotFi


    def solve(self):

        #We know E,W, N, S boundary conditions, so lets speed up things
        self.setBoundaryConditions()

        self.callTDMA()

        #print the solution vector
        print self.Fi


    #end solve

    
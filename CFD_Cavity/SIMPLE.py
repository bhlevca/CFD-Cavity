'''
Created on Apr 8, 2012

@author: bogdan
'''

import FinVol_2D_Conv_Diff
import LidCavity
import numpy
#import the TDMA module
from thomas import *

import linalg


# Create a mesh class that holds a vector of nodes


class SIMPLE(object):
    '''
    classdocs
    
    The NS equation is solved for Fi = u  and Fi = v and Gamma = miu
    Solve the equation:
    d/dx(rho*u*Fi) + d/dy(ro*v*Fi) = d/dx(Gamma* dFi/dx)+d/dy(Gamma* dFi/dy) - dp/dx +Su
    
    d/dx(rho*u*u) + d/dy(ro*v*u) = d/dx(Gamma* du/dx)+d/dy(Gamma* du/dy) - dp/dx +Su
    d/dx(rho*u*v) + d/dy(ro*v*v) = d/dx(Gamma* dv/dx)+d/dy(Gamma* dv/dy) - dp/dx +Su
    
    '''



    def __init__(self, lidCav, scheme, debug):
        '''
        Constructor
        '''
        #specific to SIMPLE
        self.itermax = 4000
        self.itmaxSIMPLE = 89 #41
        self.eps = 1e-8    #TDMA
        self.err = 1.73       #1e-3          #SIMPLE 


        #inherited from the Cavity problem
        self.lc = lidCav
        self.scheme = scheme
        self.debug = debug
        self.Nx = self.lc.Nx     #number of nodes
        self.Ny = self.lc.Ny
        self.Lx = self.lc.Lx     #domain dimension
        self.Ly = self.lc.Ly
        self.Gx = self.lc.miuX   # diffusive coeff Gamma or Miu 
        self.Gy = self.lc.miuY
        self.rho = self.lc.rho   # density
        self.deltaX = None
        self.deltaY = None
        self.finvol_u = None
        self.finvol_v = None
        self.W = self.lc.W
        self.E = self.lc.E
        self.N = self.lc.N
        self.S = self.lc.S
        self.P = self.lc.P


        if scheme == "UD":
            self.urf_p = 0.2  #underrelaxation factor for P  GOOD for 50x50=0.1   129x129=0.1     25x25=0.1 self.urf = 0.2  
            self.urf_uv = 0.6  #underrelaxation factor for U *V  GOOD for 50x50=0.6  129x129=0.6  25x25=0.6  self.urf = 0.2   
        elif scheme == 'QUICK':
            self.urf_p = 0.061 #0.061  #underrelaxation factor for P  GOOD for 25x25=0.061    
            self.urf_uv = 0.9  #underrelaxation factor for U *V  GOOD for 25x25=0.9   
        else :
            self.urf_p = 0.1  #  
            self.urf_uv = 0.6  #   

        if self.lc.CPP == True:
            self.cls = linalg.CppBlas()



        #self.urf_p = 0.3
        #self.urf_uv = 0.2
        self.urf = 0.2   #underrelaxation factor for TDMA

        #init the nodes for pressure
        self.aW = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aE = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aN = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aS = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aP = numpy.zeros((self.lc.Ny, self.lc.Nx))

        self.aWW = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aEE = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aNN = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.aSS = numpy.zeros((self.lc.Ny, self.lc.Nx))

        #init Su sources
        self.Su = numpy.zeros((self.lc.Ny, self.lc.Nx))


        #init Sp sources
        self.SpN = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.SpS = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.SpE = numpy.zeros((self.lc.Ny, self.lc.Nx))
        self.SpW = numpy.zeros((self.lc.Ny, self.lc.Nx))



        #Convection terms
        self.F = numpy.zeros((4, self.lc.Ny, self.lc.Nx))
        self.FOld = numpy.zeros((4, self.lc.Ny, self.lc.Nx))



        self.V = numpy.zeros(5)   #Volume on 4 directions + P

        self.u = numpy.zeros((self.lc.Ny, self.lc.Nx)) # horizontal velocity
        self.ue = numpy.zeros((self.lc.Ny, self.lc.Nx)) # horizontal E face velocity
        self.uw = numpy.zeros((self.lc.Ny, self.lc.Nx)) # horizontal W face velocity

        self.v = numpy.zeros((self.lc.Ny, self.lc.Nx)) # vertical velocity
        self.vn = numpy.zeros((self.lc.Ny, self.lc.Nx)) # vertical N face velocity
        self.vs = numpy.zeros((self.lc.Ny, self.lc.Nx)) # vertical S face  velocity

        self.p = numpy.zeros((self.lc.Ny, self.lc.Nx)) # pressure

        self.ustar = numpy.zeros((self.lc.Ny, self.lc.Nx)) # horizontal velocity - initial
        self.vstar = numpy.zeros((self.lc.Ny, self.lc.Nx)) # vertical velocity - initial
        self.pstar = numpy.zeros((self.lc.Ny, self.lc.Nx)) # pressure - initial

        self.pnot = numpy.zeros((self.lc.Ny, self.lc.Nx)) # pressure correction
        self.pnotOld = numpy.zeros((self.lc.Ny, self.lc.Nx)) # pressure correction

        #the TDMA coeficients. (We calculate vertical lines (columns = > allocate the number of horizontal lines)
        self.alp = numpy.zeros(self.Ny)
        self.bet = numpy.zeros(self.Ny)
        self.D = numpy.zeros(self.Ny)
        self.C = numpy.zeros(self.Ny)


        self.calculateDelta()
    #end __init__

    #calculate the sources for each face
    def Sxe(self, j, i):
        #return 0.5 * (self.finvol_u.SuP[j, i] + self.finvol_u.SuE[j, i])
        return 0.0

    def Sxw(self, j, i):
        return 0.0

    def Syn(self, j, i):
        return 0.0

    def Sys(self, j, i):
        return 0.0

    def SxE(self, j, i):
        return 0.0

    def SxW(self, j, i):
        return 0.0

    def SyN(self, j, i):
        return 0.0

    def SyS(self, j, i):
        return 0.0

    def SxP(self, j, i):
        return 0.0

    def SyP(self, j, i):
        return 0.0

    def pE(self, p, j, i):
        if i == self.Nx - 1:
            return 2 * p[j, i] - p[j, i - 1]  #extrapolate from the neighboring values
        return p[j, i + 1]

    def pW(self, p, j, i):
        if i == 0:
            return 2 * p[j, i] - p[j, i + 1]  #extrapolate from the neighboring values
        return p[j, i - 1]

    def pEE(self, p, j, i):
        if i == self.Nx - 1:
            return 0.0
        if i == self.Nx - 2:
            return 2 * p[j, i + 1] - p[j, i]  #extrapolate from the neighboring values

        return p[j, i + 2]

    def pWW(self, p, j, i):
        if i == 0:
            return 0.0
        if i == 1:
            return 2 * p[j, i - 1] - p[j, i] #extrapolate from the neighboring values
        return p[j, i - 2]

    def pN(self, p, j, i):
        if j == self.Ny - 1:
            return 2 * p[j, i] - p[j - 1, i]  #extrapolate from the neighboring values

        return p[j + 1, i]

    def pNN(self, p, j, i):
        if j == self.Ny - 1:
            return 0.0
        if j == self.Ny - 2:
            return 2 * p[j + 1, i] - p[j, i]  #extrapolate from the neighboring values

        return p[j + 2, i]

    def pS(self, p, j, i):
        if j == 0:
            return 2 * p[j, i] - p[j + 1, i]  #extrapolate from the neighboring values
        return p[j - 1, i]

    def pSS(self, p, j, i):
        if j == 0:
            return 0.0
        if j == 1:
            return 2 * p[j - 1, i] - p[j, i]  #extrapolate from the neighboring valuesself.urf_uv
        return p[j - 2, i]

    def pP(self, p, j, i):
        return p[j, i]

    def pe(self, p, j, i):
        return 0.5 * (self.pE(p, j, i) + self.pP(p, j, i))

    def pn(self, p, j, i):
        return 0.5 * (self.pN(p, j, i) + self.pP(p, j, i))

    def pw(self, p, j, i):
        return 0.5 * (self.pP(p, j, i) + self.pW(p, j, i))

    def ps(self, p, j, i):
        return 0.5 * (self.pP(p, j, i) + self.pS(p, j, i))


    def calculateDelta(self):
        '''
        calculate preliminary information, deltaX and Y 
        '''
        self.deltaX = self.Lx / self.Nx
        self.deltaY = self.Ly / self.Ny
        self.alp = numpy.zeros(self.Ny)
        self.bet = numpy.zeros(self.Ny)
        self.D = numpy.zeros(self.Ny)
        self.C = numpy.zeros(self.Ny)

        self.V[self.P] = self.deltaX * self.deltaY
        self.V[self.E] = self.deltaX * self.deltaY
        self.V[self.W] = self.deltaX * self.deltaY
        self.V[self.N] = self.deltaX * self.deltaY
        self.V[self.S] = self.deltaX * self.deltaY


    def calculate_uv_face(self):
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                self.ue[j, i] = self.ue_f(j, i)
                self.uw[j, i] = self.uw_f(j, i)
                self.vn[j, i] = self.vn_f(j, i)
                self.vs[j, i] = self.vs_f(j, i)

    def calculate_F(self, firstiter = False):
        '''
            calculate convection terms for cell faces
            unsteady flow F is NOT constant
       
        '''
        if firstiter == True:
            for i in range(0, self.Nx):
                for j in range(0, self.Ny):
                    self.F[self.E, j, i] = self.lc.rho * self.ustar[j, i] * self.V[self.E] / self.deltaY
                    self.F[self.W, j, i] = self.lc.rho * self.ustar[j, i] * self.V[self.W] / self.deltaY
                    self.F[self.N, j, i] = self.lc.rho * self.vstar[j, i] * self.V[self.N] / self.deltaX
                    self.F[self.S, j, i] = self.lc.rho * self.vstar[j, i] * self.V[self.S] / self.deltaX
        else:
            for i in range(0, self.Nx):
                for j in range(0, self.Ny):
                    self.F[self.E, j, i] = self.lc.rho * self.ue[j, i] * self.V[self.E] / self.deltaY
                    self.F[self.W, j, i] = self.lc.rho * self.uw[j, i] * self.V[self.W] / self.deltaY
                    self.F[self.N, j, i] = self.lc.rho * self.vn[j, i] * self.V[self.N] / self.deltaX
                    self.F[self.S, j, i] = self.lc.rho * self.vs[j, i] * self.V[self.S] / self.deltaX
    #end calculate_F            


    #node velocities
    #horizontal
    def uW(self, j, i):
        #FiW
        if i == 0:
            u = self.lc.FiU[self.W]
        else:
            u = self.ustar[j, i - 1]
        return u

    def uE(self, j, i):
        if i == self.Nx - 1:
            u = self.lc.FiU[self.E]
        else:
            u = self.ustar[j, i + 1]
        return u

    def uS(self, j, i):
        if j == 0:
            u = self.lc.FiU[self.S]
        else:
            u = self.ustar[j - 1, i]
        return u

    def uN(self, j, i):
        #FiW
        if j == self.Ny - 1:
            u = self.lc.FiU[self.N]
        else:
            u = self.ustar[j + 1, i]
        return u
    #vertical
    def vW(self, j, i):
        #FiW
        if i == 0:
            v = self.lc.FiV[self.W]
        else:
            v = self.vstar[j, i - 1]
        return v

    def vE(self, j, i):
        if i == self.Nx - 1:
            v = self.lc.FiV[self.E]
        else:
            v = self.vstar[j, i + 1]
        return v

    def vS(self, j, i):
        if j == 0:
            v = self.lc.FiV[self.S]
        else:
            v = self.vstar[j - 1, i]
        return v

    def vN(self, j, i):
        #FiW
        if j == self.Ny - 1:
            v = self.lc.FiV[self.N]
        else:
            v = self.vstar[j + 1, i]
        return v

    #face velocities
    #horizontal  
    def ue_f(self, j, i):
        if  self.finvol_u.aE[j, i] == 0 or i == self.Nx - 1:
            #VEaE = 0
            return self.lc.ue_wall[j]
        else:
            VEaE = self.V[self.E] / self.finvol_u.aE[j, i]
        VPaP = self.V[self.P] / self.finvol_u.aP[j, i]
        Veae = 0.5 * (VEaE + VPaP)
        #Veae = (self.V[self.E] + self.V[self.P]) / (self.finvol_u.aE[j, i] + self.finvol_u.aP[j, i])
        dpdxE = (self.pEE(self.pstar, j, i) - self.pP(self.pstar, j, i)) / (2.0 * self.deltaX)
        dpdxP = (self.pE(self.pstar, j, i) - self.pW(self.pstar, j, i)) / (2.0 * self.deltaX)
        if i == self.Nx - 1:
            dpdxe = dpdxE
        else:
            dpdxe = 0.5 * (dpdxE + dpdxP) # 
            #dpdxe = (self.pE(self.pstar, j, i) - self.pP(self.pstar, j, i)) / self.deltaX
        ue = (self.uE(j, i) + self.ustar[j, i]) / 2.0 + \
            self.urf_uv * 0.5 * (VEaE * (dpdxE - self.SxE(j, i)) + VPaP * (dpdxP - self.SxP(j, i))) - Veae * (dpdxe - self.Sxe(j, i))

        #altenative formula
        #ue2 = (self.uE(j, i) + self.ustar[j, i]) / 2.0 + \
        #    0.5 * (self.deltaY / self.finvol_u.aE[j, i] * (self.pEE(self.pstar, j, i) - self.pP(self.pstar, j, i)) / 2.0\
        #            + self.deltaY / self.finvol_u.aP[j, i] * (self.pE(self.pstar, j, i) - self.pW(self.pstar, j, i)) / 2.0)\
        #            - self.deltaY * 0.5 * (1.0 / self.finvol_u.aP[j, i] + 1.0 / self.finvol_u.aE[j, i]) * (self.pE(self.pstar, j, i) - self.pP(self.pstar, j, i))


        return ue

    def uw_f(self, j, i):
        if  self.finvol_u.aW[j, i] == 0 or i == 0 :
            #VWaW = 0
            return self.lc.uw_wall[j]
        else:
            VWaW = self.V[self.W] / self.finvol_u.aW[j, i]
        VPaP = self.V[self.P] / self.finvol_u.aP[j, i]
        Vwaw = 0.5 * (VWaW + VPaP)
        #Vwaw = (self.V[self.W] + self.V[self.P]) / (self.finvol_u.aW[j, i] + self.finvol_u.aP[j, i])
        dpdxW = (self.pP(self.pstar, j, i) - self.pWW(self.pstar, j, i)) / (2 * self.deltaX)
        dpdxP = (self.pE(self.pstar, j, i) - self.pW(self.pstar, j, i)) / (2 * self.deltaX)
        if i == 0:
            dpdxw = dpdxW
        else:
            dpdxw = 0.5 * (dpdxW + dpdxP)
            #dpdxw = (self.pP(self.pstar, j, i) - self.pW(self.pstar, j, i)) / self.deltaX
        uw = (self.uW(j, i) + self.ustar[j, i]) / 2 + \
            self.urf_uv * 0.5 * (VWaW * (dpdxW - self.SxW(j, i)) + VPaP * (dpdxP - self.SxP(j, i))) - Vwaw * (dpdxw - self.Sxw(j, i))
        #or
        #uv2 = (self.uW(j, i) + self.ustar[j, i]) / 2 + Vwaw * (0.5 * (dpdxW + dpdxP) - dpdxw)
        return uw

    #vertical
    def vn_f(self, j, i):
        if  self.finvol_u.aN[j, i] == 0 or j == self.Ny - 1 :
            #VNaN = 0
            return self.lc.vn_wall[i]
        else:
            VNaN = self.V[self.N] / self.finvol_v.aN[j, i]
        VPaP = self.V[self.P] / self.finvol_v.aP[j, i]
        Vnan = 0.5 * (VNaN + VPaP)
        #Vnan = (self.V[self.N] + self.V[self.P]) / (self.finvol_u.aN[j, i] + self.finvol_u.aP[j, i])
        dpdxN = (self.pNN(self.pstar, j, i) - self.pP(self.pstar, j, i)) / (2 * self.deltaY)
        dpdxP = (self.pN(self.pstar, j, i) - self.pS(self.pstar, j, i)) / (2 * self.deltaY)
        if j == self.Ny - 1:
            dpdxn = dpdxN
        else:
            dpdxn = 0.5 * (dpdxN + dpdxP)
            #dpdxn = (self.pN(self.pstar, j, i) - self.pP(self.pstar, j, i)) / self.deltaY
        vn = (self.vN(j, i) + self.vstar[j, i]) / 2 + \
            self.urf_uv * 0.5 * (VNaN * (dpdxN - self.SyN(j, i)) + VPaP * (dpdxP - self.SxP(j, i))) - Vnan * (dpdxn - self.Syn(j, i))
        return vn


    def vs_f(self, j, i):
        if  self.finvol_u.aS[j, i] == 0 or j == 0:
            #VSaS = 0
            return self.lc.vs_wall[i]
        else:
            VSaS = self.V[self.S] / self.finvol_v.aS[j, i]
        VPaP = self.V[self.P] / self.finvol_v.aP[j, i]
        Vsas = 0.5 * (VSaS + VPaP)
        #Vsas = (self.V[self.S] + self.V[self.P]) / (self.finvol_u.aS[j, i] + self.finvol_u.aP[j, i])
        dpdxS = (self.pP(self.pstar, j, i) - self.pSS(self.pstar, j, i)) / (2 * self.deltaY)
        dpdxP = (self.pN(self.pstar, j, i) - self.pS(self.pstar, j, i)) / (2 * self.deltaY)
        if j == 0:
            dpdxs = dpdxS
        else:
            dpdxs = 0.5 * (dpdxS + dpdxP)
            #dpdxs = (self.pP(self.pstar, j, i) - self.pS(self.pstar, j, i)) / self.deltaY
        vs = (self.vS(j, i) + self.vstar[j, i]) / 2 + \
            self.urf_uv * 0.5 * (VSaS * (dpdxS - self.SyS(j, i)) + VPaP * (dpdxP - self.SxP(j, i))) - Vsas * (dpdxs - self.Sys(j, i))
        return vs



    def solve(self):
        converged = False
        Iter = 0

        # 1) Estimate intial values for the ustar, vstar and pstar variables
        self.guess_p_u_v()

        # Calcuate intitial convection terms
        self.calculate_F(Iter == 0)

        while  converged == False and Iter < self.itmaxSIMPLE:
            # 2) solve the momentum equations for the star/initial values
            self.finvol_u, self.finvol_v = self.solve_momentum_u_v()

            # 3) recalculate interface velocity and new convective flux terms (Fs) 
            #if Iter == 0 :
            self.calculate_uv_face()
            self.calculate_F()

            # 4) solve aP * p'P = aW * p'W + aE * p'E  + aS*p'S + aN* p'N + b
            self.solve_pnot(self.finvol_u, self.finvol_v)


            # 8) Extrapolate to boundaries
            self.extrapolate_p_to_boundaries()

            # 5) correct 
            #     pressure p = p* + p'
            self.correct_p_from_pnot_pstar()
            #     corect velocities u = u* + (dP*u')
            self.correct_uv()

            # 6) correct m not - Convection terms ex: Fe= Fe +( rho *de * deltaY)*(Pp' -Pe')
            # => this step seems to make convergence worse
            #self.correct_convection_terms()

            # 7) optional - solve other transports
            #self.solve_other_transports()



            # 9) check convergence
            converged = self.check_convergence(Iter)

            # 10) set the new values 
            self.set_new_values()

            # inrement the Iterator
            Iter += 1

        if Iter > self.itermax:
            print "Solution did not converge in %d interations", Iter
        else:
            print "Solution converged in %d interations", Iter
        return [self.u, self.v]
    #end solve

    def guess_p_u_v(self):
        '''
            Estimate intial valuess for the u, v and p variables
        '''

        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                self.ustar[j, i] = 0
                self.vstar[j, i] = 0
                self.pstar[j, i] = 0 #9.8 * (self.Ny - 1 - j) * self.deltaY * self.rho
                #other may follow , Temp, Conc, etc
    #end guess_p_u_v


    def solve_momentum_u_v(self):
        #set boundary conditions for Fi = v , all are 0 
        self.finvol_v = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff("v", self.lc.FiV, self.pstar, self)
        self.vstar = self.finvol_v.solve().copy()


        #set boundary conditions for Fi = u 
        self.finvol_u = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff("u", self.lc.FiU, self.pstar, self)
        self.ustar = self.finvol_u.solve().copy()
        return [self.finvol_u, self.finvol_v]



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
            self.C[j] = self.aW[j, i] * self.pW(self.pnot, j, i) + self.aE[j, i] * self.pE(self.pnot, j, i) + self.Su[j, i]
        #end for j self.solve_pnot(finvol_u, finvol_v)



    #end calculateTDMACoefficients


    def solve_pnot(self, finvol_u, finvol_v):
        #calculate the coefficients for pnot 
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):

                if i == self.Nx - 1:
                    self.aE[j, i] = 0
                else:
                    ape = 0.5 * (self.finvol_u.aE[j, i] + self.finvol_u.aP[j, i])
                    #de = self.finvol_u.A[self.E] / ape
                    de = 1.0 / ape
                    #de = 0.5 * (self.finvol_u.A[self.E] / self.finvol_u.aE[j, i] + self.finvol_u.A[self.E] / self.finvol_u.aP[j, i])
                    self.aE[j, i] = self.rho * de * self.deltaY

                if i == 0:
                    self.aW[j, i] = 0
                else:
                    apw = 0.5 * (self.finvol_u.aW[j, i] + self.finvol_u.aP[j, i])
                    #dw = self.finvol_u.A[self.W] / apw
                    dw = 1.0 / apw
                    #dw = 0.5 * (self.finvol_u.A[self.W] / self.finvol_u.aW[j, i] + self.finvol_u.A[self.W] / self.finvol_u.aP[j, i])
                    self.aW[j, i] = self.rho * dw * self.deltaY

                if j == self.Ny - 1:
                    self.aN[j, i] = 0
                else:
                    apn = 0.5 * (self.finvol_v.aN[j, i] + self.finvol_v.aP[j, i])
                    #dn = self.finvol_v.A[self.N] / apn
                    dn = 1.0 / apn
                    #dn = 0.5 * (self.finvol_v.A[self.N] / self.finvol_v.aN[j, i] + self.finvol_v.A[self.N] / self.finvol_v.aP[j, i])
                    self.aN[j, i] = self.rho * dn * self.deltaX

                if j == 0:
                    self.aS[j, i] = 0
                else:
                    aps = 0.5 * (self.finvol_v.aS[j, i] + self.finvol_v.aP[j, i])
                    #ds = self.finvol_v.A[self.S] / aps
                    ds = 1.0 / aps
                    #ds = 0.5 * (self.finvol_v.A[self.E] / self.finvol_v.aS[j, i] + self.finvol_v.A[self.S] / self.finvol_v.aP[j, i])
                    self.aS[j, i] = self.rho * ds * self.deltaX

                self.Su[j, i] = (self.rho * self.uw[j, i] - self.rho * self.ue[j, i]) * self.deltaY + \
                                (self.rho * self.vs[j, i] - self.rho * self.vn[j, i]) * self.deltaX

                self.aP[j, i] = self.aE[j, i] + self.aW[j, i] + self.aN[j, i] + self.aS[j, i]

        it = 0
        n = self.D.size
        x = numpy.zeros(n)

        while  self.itermax > it :
            #copy current values to the old values matrix
            self.pnotOld = self.pnot.copy()

            #Swipe  from W to E
            for i in range(0, self.Nx):

                #calculate the TDMA coefficients for column i
                self.calculateTDMACoefficients(i)

                if self.debug == True:
                    print "beta:", self.bet
                    print "D", self.D
                    print "alp", self.alp
                    print "C", self.C

                if self.lc.CPP == True:
                    self.cls.setTDMA(-self.bet[1:], self.D, -self.alp[:-1], self.C, n)
                    d = self.cls.solveTDMA(x, n)
                    self.pnot[:, i] = d["solution"].copy()
                else:
                    x = thomas(n, -self.bet[1:], self.D, -self.alp[:-1], self.C)
                    self.pnot[:, i] = x.copy()

            #end i

            #TODO Under relaxation
            self.pnot = self.urf * self.pnot.copy() + self.pnotOld.copy() * (1 - self.urf)

            #test accuracy and exit condition
            flat = self.pnot[:, 1] - self.pnotOld[:, 1]
            dx = math.sqrt(numpy.dot(flat, flat))

            if it % 600 == 0:
                print "var: %s iter # %d, dx=%1.9f" % ("p", it, dx)
            #print "Fi:", self.Fi

            #Exit if we are satisfied wit the accuracy
            if dx < self.eps :
                print self.pnot
                return

            it += 1

        #end while

        #if we did not converge yet print an error and exit
        if self.itermax <= it:
            print "Max iterations exceeded => did not converge"
            print self.pnot
            return


        return [self.pnot]
    #end solve_pnot


    def correct_p_from_pnot_pstar(self):
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                self.p[j, i] = self.pstar[j, i].copy() + self.urf_p * self.pnot[j, i]

    def correct_uv(self):
        '''
        correct u /v from u*/v*  +  p' 
        '''
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                self.u[j, i] = self.ustar[j, i].copy() + self.urf_uv * self.finvol_u.A[self.E] / self.finvol_u.aP[j, i] * (self.pw(self.pnot, j, i) - self.pe(self.pnot, j, i))
                self.v[j, i] = self.vstar[j, i].copy() + self.urf_uv * self.finvol_v.A[self.N] / self.finvol_v.aP[j, i] * (self.ps(self.pnot, j, i) - self.pn(self.pnot, j, i))
    #end correct_uv

    def update_F(self):
        self.FOld = self.F.copy()
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                self.F[self.E, j, i] = self.FOld[self.E, j, i].copy() + self.rho * self.ue[j, i] * self.deltaY
                self.F[self.W, j, i] = self.FOld[self.W, j, i].copy() + self.rho * self.uw[j, i] * self.deltaY
                self.F[self.N, j, i] = self.FOld[self.N, j, i].copy() + self.rho * self.vn[j, i] * self.deltaX
                self.F[self.S, j, i] = self.FOld[self.S, j, i].copy() + self.rho * self.vs[j, i] * self.deltaX

    def correct_convection_terms(self):
        self.calculate_uv_face()
        self.calculate_F()
        return
        for i in range(0, self.Nx):
            for j in range(0, self.Ny):
                if i == self.Nx - 1:
                    self.ue[j, i] = 0
                else:
                    de = self.urf_uv * 0.5 * (self.finvol_u.A[self.E] / self.finvol_u.aE[j, i] + self.finvol_u.A[self.E] / self.finvol_u.aP[j, i])
                    #de = self.urf_uv * 0.5* (1.0 / self.finvol_u.aE[j, i] + 1.0 / self.finvol_u.aP[j, i])
                    self.ue[j, i] = self.ue[j, i] + de * (self.pP(self.pnot, j, i) - self.pE(self.pnot, j, i))

                if i == 0:
                    self.uw[j, i] = 0
                else:
                    dw = self.urf_uv * 0.5 * (self.finvol_u.A[self.W] / self.finvol_u.aW[j, i] + self.finvol_u.A[self.W] / self.finvol_u.aP[j, i])
                    #dw = self.urf_uv * 0.5* (1.0 / self.finvol_u.aW[j, i] + 1.0 / self.finvol_u.aP[j, i])
                    self.uw[j, i] = self.uw[j, i] + dw * (self.pW(self.pnot, j, i) - self.pP(self.pnot, j, i))

                if j == self.Ny - 1:
                    self.vn[j, i] = 0
                else:
                    dn = self.urf_uv * 0.5 * (self.finvol_v.A[self.N] / self.finvol_v.aN[j, i] + self.finvol_v.A[self.W] / self.finvol_u.aP[j, i])
                    #dn = self.urf_uv * 0.5* (1.0 / self.finvol_u.aN[j, i] + 1.0 / self.finvol_u.aP[j, i])
                    self.vn[j, i] = self.vn[j, i] + dn * (self.pP(self.pnot, j, i) - self.pN(self.pnot, j, i))

                if j == 0:
                    self.vs[j, i] = 0
                else:
                    ds = self.urf_uv * 0.5 * (self.finvol_v.A[self.S] / self.finvol_v.aS[j, i] + self.finvol_v.A[self.W] / self.finvol_u.aP[j, i])
                    #ds = self.urf_uv * 0.5* (1.0 / self.finvol_u.aE[j, i] + 1.0 / self.finvol_u.aP[j, i])
                    self.vs[j, i] = self.vs[j, i] + ds * (self.pS(self.pnot, j, i) - self.pP(self.pnot, j, i))

        self.update_F()

    def solve_other_transports(self):
        '''
            Could be temperature, Concentration, etc.
            Do nothing for now.
        '''
        pass

    def extrapolate_p_to_boundaries(self):
        '''
            The pressure values at the boundary conditions can be calculated by linear
            interpolation of using the the two near boundary node pressures
        '''
        # this is done automatically in the pE functions above for p but not for p'
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                if i == self.Nx - 1:
                    self.p[j, i] = self.p[j, i - 1]
                    self.pnot[j, i] = self.pnot[j, i - 1];
                if i == 0:
                    self.p[j, i] = self.p[j, i + 1]
                    self.pnot[j, i] = self.pnot[j, i + 1];
                if j == self.Ny - 1:
                    self.p[j, i] = self.p[j - 1, i]
                    self.pnot[j, i] = self.pnot[j - 1, i];
                if j == 0:
                    self.p[j, i] = self.p[j + 1, i]
                    self.pnot[j, i] = self.pnot[j + 1, i];

    #end extrapolate_p_to_boundaries

    def check_convergence(self, Iter):
        '''
            All 3 variables need to converge
            
            
        '''
        Sum = 0
        for j in range(0, self.Ny):
            for i in range(0, self.Nx):
                Sum += abs(self.F[self.E, j, i] - self.F[self.W, j, i] + self.F[self.N, j, i] - self.F[self.S, j, i])

        print "iteration # %d, Sum abs(F) =%1.9f  " % (Iter, Sum)

        #Exit if we are satisfied wit the accuracy
        if Sum <= self.err :
            return True
        return False

    def set_new_values(self):
        self.pstar = self.p.copy()
        self.ustar = self.u.copy()
        self.vstar = self.v.copy()







'''
Created on Feb 7, 2012

@author: bogdan hlevca

Diffusion on a square domain , Written in Python and numerical Python using matplotlib library
'''

#import necessary libraries
import numpy

#import the program
import FinVol_2D_Diffusion
import FinVol_2D_DiffusionEx7_2


#start main program in Python
if __name__ == '__main__':

    '''
    Solve the equation:
    d/dx(Gamma*A* dF/dx)+d/dy(Gamma*A* dF/dy) + h*L(F-F0) = 0 
    
    where 
        d is the partial difference d_rond
        F - stands for the property Fi
    '''

    #Given values 
    Lx = 1.0            #domain rectangle length on X axis [m]
    Ly = 1.0            #domain rectangle length on Y axis                 
    Gx = 20.0           # Diffusion coefficient on X axis  [W/(m^2 * K)]
    Gy = 20.0           # Diffusion coefficient on X axis  
    h = 10.0            # convective heat transfer coefficient [W/(m^2 * K)] 
    Fi0 = 300.0         # property at the in the N boundary flux [C]
    FiE = 100.0         # property at the E boundary [C]
    FiW = 10.0          # property at the W boundary [C] 
    Nx = 80             # Number of nodes on X axis   [10 , 20, 40 ,80]
    Ny = 80             # Number of nodes on Y axis   [10 , 20, 40 ,80]
    debug = False       # debug flag controls the printout verbosity


    #initialize the problem
    finvol = FinVol_2D_Diffusion.FinVol_2D_Diffusion(Lx, Ly, Gx, Gy, h, Fi0, FiE, FiW, Nx, Ny, debug)

    #solve the problem
    finvol.solve()

    #plot the results
    finvol.plotFi()


    #set to True to Solve the problem 7.2 from Versteeg & Malalssekera, Introduction to CFD:the finite volume method, Second Edition(2007) 
    if False:
        Lx = 0.3
        Ly = 0.4
        Gx = 1000
        Gy = 1000

        FiN = 100.0
        FiE = 100.0
        FiW = 10.0
        Nx = 3
        Ny = 4
        qw = 500 * 1000
        debug = True
        d = 0.01


        finvol = FinVol_2D_DiffusionEx7_2.FinVol_2D_Diffusion_Ex7_2(Lx, Ly, Gx, Gy, FiN, Nx, Ny, qw, d, True)
        result = finvol.solveEx()


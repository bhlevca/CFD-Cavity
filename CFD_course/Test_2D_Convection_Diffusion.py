'''
Created on March 7, 2012

@author: bogdan hlevca 995151213

Diffusion on a square domain , Written in Python and numerical Python using matplotlib library
'''

#import necessary libraries
import numpy

#import the program
import FinVol_TVD_2D_Conv_Diff
import FinVol_2D_Conv_Diff
import FinVol_QUICK_Hayase_2D_Conv_Diff


def diag2(arr):
    n = arr.shape[0]
    d = numpy.zeros(n)
    for j in range(0, n):
        for i in range(0, n):
            if j == n - 1 - i:
                d[j] = arr[j, i]
        #end for
    #end for
    reversed_arr = d[::-1]
    return reversed_arr

#start main program in Python
if __name__ == '__main__':

    '''
    Solve the equation:
    d/dx(Gamma*A* dF/dx)+d/dy(Gamma*A* dF/dy) + h*L(F-F0) = 0 
    
    where 
        d is the partial difference d_rond
        F - stands for the property Fi
    '''


    test = 'Compare' #False
    test = False

    #Given values
    if test == True:
        #This is the 5.5 example from Versteeg 2007 using QUICK method
        scheme = 'QUICK'
        Lx = 1.0            #domain rectangle length on X axis [m]
        Ly = 1.0            #domain rectangle length on Y axis                 
        Gx = 0.1           # Diffusion coefficient on X axis  [W/(m^2 * K)]
        Gy = 0.1           # Diffusion coefficient on X axis  
        FiE = 0.0         # property at the E boundary [C]
        FiW = 0.0          # property at the W boundary [C] 
        FiN = 0.0         # property at the E boundary [C]
        FiS = 1.0          # property at the W boundary [C] 
        rho = 1.0        # density 100 Kg/m^3
        Nx = 1             # Number of nodes on X axis   [10 , 40 , 100, 300]
        Ny = 5           # Number of nodes on Y axis   [10 , 40,  100 , 300]
        u = 0.2              # W->E velocity 2m/s 
        v = 0.2             # S->N velocity 2m/s
        debug = True       # debug flag controls the printout verbosity
    else:
        Lx = 1.0            #domain rectangle length on X axis [m]
        Ly = 1.0            #domain rectangle length on Y axis                 
        Gx = 1.0           # Diffusion coefficient on X axis  [W/(m^2 * K)]
        Gy = 1.0           # Diffusion coefficient on X axis  
        FiE = 0.0         # property at the E boundary [C]
        FiW = 100.0          # property at the W boundary [C] 
        FiN = 100.0         # property at the E boundary [C]
        FiS = 0.0          # property at the W boundary [C] 
        rho = 100.0        # density 100 Kg/m^3
        Nx = 40             # Number of nodes on X axis   [10 , 40 , 100, 300]
        Ny = 40           # Number of nodes on Y axis   [10 , 40,  100 , 300]
        u = 2.0              # W->E velocity 2m/s 
        v = 2.0             # S->N velocity 2m/s
        debug = False       # debug flag controls the printout verbosity
        scheme = 'QUICK'
        #scheme = 'CD'
        #scheme = 'UD'
    #endif



    #initialize the problem

    if test == True:
        #finvol = FinVol_TVD_2D_Conv_Diff.FinVol_TVD_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, debug)
        finvol = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, scheme, debug)
        #solve
        finvol.solve()
        #plot the results
        finvol.plotFi()

    elif test == 'Compare':
        finvol1 = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, "CD", debug)
        FiCD = finvol1.solve()
        finvol2 = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, "UD", debug)
        FiUD = finvol2.solve()
        finvol3 = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, "QUICK", debug)
        FiQU = finvol3.solve()
        #plot the 3 diagonals



        finvol3.plot3curves(diag2(FiCD), diag2(FiUD), diag2(FiQU))
    else:
        #finvol = FinVol_QUICK_Hayase_2D_Conv_Diff.FinVol_QUICK_Hayase_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, debug)
        finvol = FinVol_2D_Conv_Diff.FinVol_2D_Conv_Diff(Lx, Ly, Gx, Gy, FiE, FiW, FiN, FiS, Nx, Ny, u, v, rho, scheme, debug)
        #solve
        finvol.solve()
        #plot the results
        finvol.plotFi()

    #solve the problem





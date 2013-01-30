'''
Created on March 7, 2012

@author: bogdan hlevca 995151213

Diffusion on a square domain , Written in Python and numerical Python using matplotlib library
'''

#import necessary libraries
import numpy

#import the program
import FinVol_2D_Conv_Diff
import LidCavity
import time

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
    Lid driven Cavity classical CFD testing problem
    '''

    debug = False       # debug flag controls the printout verbosity
    scheme = 'QUICK'
    #scheme = 'CD'
    scheme = 'UD'

    cavity = LidCavity.LidCavity(scheme, debug)

    start = time.time()

    u, v = cavity.solve()
    #plot the results
    cavity.plotStreams(u, v)
    cavity.plotCentreLinesProfiles(u, v)

    end = time.time()
    elapsed = end - start
    print "Time taken: ", elapsed, "seconds."

    #cavity.plotHeatMap(u, v)
    cavity.showGraph()






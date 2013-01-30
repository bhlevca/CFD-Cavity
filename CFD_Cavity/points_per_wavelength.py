import numpy
import math

#import the graphi library
import matplotlib.pyplot as plt

eps = 0.001

kh = numpy.arange(0., 3.0, 0.1)              #create the x, and y divisions
kstarh = numpy.zeros((kh.size))

for i in range(0, kh.size):
    kstarh[i] = (4 / 3 * math.sin(kh[i]) - 1 / 6 * math.sin(2 * kh[i]))


#draw the 
fig2 = plt.figure(2, facecolor = 'w', edgecolor = 'k')   #prepare the plotting environment
ax = fig2.add_subplot(111)
im = ax.plot(kh, kh, kh, kstarh)
plt.show()


#find the value
for i in range(0, kh.size):
    if i > 0 and abs((kh[i] - kstarh[i]) / kh[i]) > eps:
        sol = kh[i + 1]
        break



ppw = 2 * math.pi / sol / 2

print "ppw = %d" % ppw

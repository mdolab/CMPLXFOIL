#! /usr/bin/env python
import pyxlight as xfoil
import pyxlight_cs as xfoil_cs
import string
from numpy import zeros

# Opterating Condition
RE = 1.0e6 # Reynold number = 1,000,000
M = .02    # Mach Number = 0.02
ITMAX = 50 # Maximum of 50 Iterations

# Load a naca0012 airfoil. There are 160 points in file
f = open('naca0012.dat','r')
NB = 160
x = zeros(NB)
y = zeros(NB)

i=0
for line in f:
    x[i] = float(string.split(line)[0])
    y[i] = float(string.split(line)[1])
    i = i + 1

N = 572 # This is VERY Important: The airfoil input must be a FIXED
         # length of 572. Simple set the coordinates up to NB and
         # leave the remainder as zeros

x_input = zeros([N])
y_input = zeros([N])

x_input[0:NB] = x
y_input[0:NB] = y

#Set real XFOIL data

xfoil.ci04.nb = NB       # Number of points in airfoil we are using
xfoil.cr14.xb = x_input  # X-coordinates of airfoil points we are using
xfoil.cr14.yb = y_input  # Y-coordinates of airfoil points we are using
xfoil.cr15.reinf1 = RE   # Reynolds number 
xfoil.cr09.minf1 = M     # Mach Number set
xfoil.ci04.itmax = ITMAX # Iterations Limit Set

xfoil.xfoil()       # This function is called to 'set' the airfoil buffer 

#Set complex XFOIL data

xfoil_cs.ci04.nb = NB       # Number of points in airfoil we are using
xfoil_cs.cr14.xb = x_input  # X-coordinates of airfoil points we are using
xfoil_cs.cr14.yb = y_input  # Y-coordinates of airfoil points we are using
xfoil_cs.cr15.reinf1 = RE   # Reynolds number 
xfoil_cs.cr09.minf1 = M     # Mach Number set
xfoil_cs.ci04.itmax = ITMAX # Iterations Limit Set

xfoil_cs.xfoil()       # This function is called to 'set' the airfoil buffer 

# Do an alpha sweep 
for i in xrange(10):
    xfoil.cr09.adeg = i/1.0
    xfoil.oper()
    print 'a,cl,cd (real):',i,xfoil.cr09.cl, xfoil.cr09.cd # Extract lift and drag

    xfoil_cs.cr09.adeg = i/1.0
    xfoil_cs.oper()
    print 'a,cl,cd (complex) :',i,xfoil_cs.cr09.cl, xfoil_cs.cr09.cd # Extract lift and drag


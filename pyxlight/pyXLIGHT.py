"""
pyXLIGHT

pyXLIGHT is a wrapper for Mark Drela's Xfoil code. The purpose of this
class is to provide an easy to use wrapper for xfoil for intergration
into other projects. Both real and complex versions of xfoil can be used.

Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2009)
"""

__version__ = "$Revision: $"

# =============================================================================
# Standard Python modules
# =============================================================================

import os, sys, copy, pdb, time

# =============================================================================
# External Python modules
# =============================================================================

from numpy import zeros, array

# =============================================================================
# Extension modules
# =============================================================================

from . import libxlight as xfoil
from . import libxlight_cs as xfoil_cs


class xfoilAnalysis:
    def __init__(self, airfoil_file, re=1e5, mach=0.0, iter=50):
        try:
            f = open(airfoil_file, "r")
        except:
            print("There was an error opening the airfoil file %s" % (airfoil_file))
            sys.exit(1)
        # end if

        self.re = re
        self.mach = mach
        self.iter = iter

        # Read the airfoil file
        x = []
        y = []
        for line in f:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))
        # end if
        self.x = array(x)
        self.y = array(y)

        self.setCoordinates(self.x, self.y)
        self.setCoordinatesComplex(self.x, self.y)

        return

    def setCoordinates(self, x, y):
        """Set the coordinates x,y into xfoil"""

        N = 572  # This is VERY Important: The airfoil input must be a
        # FIXED length of 572. Simply set the coordinates up
        # to NB and leave the remainder as zeros
        NB = len(x)
        x_input = zeros(N)
        y_input = zeros(N)
        x_input[:NB] = array(x).copy()
        y_input[:NB] = array(y).copy()
        xfoil.ci04.nb = NB
        xfoil.cr14.xb = x_input
        xfoil.cr14.yb = y_input
        xfoil.xfoil()

        return

    def setCoordinatesComplex(self, x, y):

        N = 572  # This is VERY Important: The airfoil input must be a
        # FIXED length of 572. Simply set the coordinates up
        # to NB and leave the remainder as zeros
        NB = len(x)
        x_input = zeros(N, "D")
        y_input = zeros(N, "D")
        x_input[:NB] = array(x).copy()
        y_input[:NB] = array(y).copy()
        xfoil_cs.ci04.nb = NB
        xfoil_cs.cr14.xb = x_input
        xfoil_cs.cr14.yb = y_input
        xfoil_cs.xfoil()

        return

    def solveAlpha(self, angle):
        """Compute the flow solution at an angle of attack angle (degrees)"""

        xfoil.cr15.reinf1 = self.re  # Reynolds number
        xfoil.cr09.minf1 = self.mach  # Mach Number set
        xfoil.ci04.itmax = self.iter  # Iterations Limit Set
        xfoil.cr09.adeg = angle
        xfoil.oper()
        return xfoil.cr09.cl, xfoil.cr09.cd, xfoil.cr09.cm, xfoil.cl01.lexitflag

    def solveAlphaComplex(self, angle):
        """Compute the flow solution at and angle of attack angle"""

        xfoil_cs.cr15.reinf1 = self.re  # Reynolds number
        xfoil_cs.cr09.minf1 = self.mach  # Mach Number set
        xfoil_cs.ci04.itmax = self.iter  # Iterations Limit Set
        xfoil_cs.cr09.adeg = angle
        xfoil_cs.oper()
        return xfoil_cs.cr09.cl, xfoil_cs.cr09.cd, xfoil_cs.cr09.cm, xfoil_cs.cl01.lexitflag

    def setValue(self, common_block, variable, value):
        """
        Set any xfoil (real) option according to:
        xfoil.<common_block>.<variable> = <value>
        """
        # First check that the attribute the user is trying to set
        # actually exists. These modules work such that you can
        # dynamically create new attributes thus if the varibale
        # doesn't exist it is just created...which is not what we want
        exec("has_attribute = hasattr(xfoil.%s,'%s')" % (common_block, variable))
        if has_attribute:
            # Now we can set it
            exec("xfoil.%s.%s = value" % (common_block, variable))
        else:
            print("There was an error in setValue")
            print("Either %s or %s does not exist" % (common_block, variable))
        # end if
        return

    def setValueComplex(self, common_block, variable, value):
        """
        Set any xfoil (complex) option according to:
        xfoil.<common_block>.<variable> = <value>
        """
        # First check that the attribute the user is trying to set
        # actually exists. These modules work such that you can
        # dynamically create new attributes thus if the varibale
        # doesn't exist it is just created...which is not what we want
        exec("has_attribute = hasattr(xfoil_cs.%s,'%s')" % (common_block, variable))
        if has_attribute:
            # Now we can set it
            exec("xfoil_cs.%s.%s = value" % (common_block, variable))
        else:
            print("There was an error in setValueComplex")
            print("Either %s or %s does not exist" % (common_block, variable))
        # end if
        return

    def xdriver(self):
        cl, cd = xfoil.xdriver(self.x, self.y)
        print("Cl=", cl, "Cd=", cd)

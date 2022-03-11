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

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
import os
import time

# =============================================================================
# Extension modules
# =============================================================================
from . import MExt


class xfoilAnalysis:
    """
    Create an xfoilAnalysis object.

    Parameters
    ----------
    airfoil_file : str
        File with list of airfoil coordinates (in .dat file format)
    re : float, optional
        Reynolds number, by default 1e5
    mach : float, optional
        Mach number, by default 0.0
    iter : int, optional
        Maximum iterations for XFOIL solver, by default 50
    debug : bool, optional
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging, by default False
    """
    def __init__(self, airfoil_file, re=1e5, mach=0.0, iter=50, debug=False):  # noqa: A002
        try:
            f = open(airfoil_file, "r")
        except OSError:
            raise OSError("There was an error opening the airfoil file %s" % (airfoil_file))
        # end if

        # Load the compiled module using MExt, allowing multiple imports
        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        time.sleep(.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil = MExt.MExt("libxlight", curDir, debug=debug)._module
        time.sleep(.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil_cs = MExt.MExt("libxlight_cs", curDir, debug=debug)._module

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
        self.x = np.array(x)
        self.y = np.array(y)

        self.setCoordinates(self.x, self.y)
        self.setCoordinatesComplex(self.x, self.y)

        return

    def setCoordinates(self, x, y):
        """Set the coordinates x,y into xfoil"""

        N = 572  # This is VERY Important: The airfoil input must be a
        # FIXED length of 572. Simply set the coordinates up
        # to NB and leave the remainder as zeros
        NB = len(x)
        x_input = np.zeros(N)
        y_input = np.zeros(N)
        x_input[:NB] = np.array(x).copy()
        y_input[:NB] = np.array(y).copy()
        self.xfoil.ci04.nb = NB
        self.xfoil.cr14.xb = x_input
        self.xfoil.cr14.yb = y_input
        self.xfoil.xfoil()

        return

    def setCoordinatesComplex(self, x, y):

        N = 572  # This is VERY Important: The airfoil input must be a
        # FIXED length of 572. Simply set the coordinates up
        # to NB and leave the remainder as zeros
        NB = len(x)
        x_input = np.zeros(N, "D")
        y_input = np.zeros(N, "D")
        x_input[:NB] = np.array(x).copy()
        y_input[:NB] = np.array(y).copy()
        self.xfoil_cs.ci04.nb = NB
        self.xfoil_cs.cr14.xb = x_input
        self.xfoil_cs.cr14.yb = y_input
        self.xfoil_cs.xfoil()

        return

    def solveAlpha(self, angle):
        """Compute the flow solution at an angle of attack angle (degrees)"""

        self.xfoil.cr15.reinf1 = self.re  # Reynolds number
        self.xfoil.cr09.minf1 = self.mach  # Mach Number set
        self.xfoil.ci04.itmax = self.iter  # Iterations Limit Set
        self.xfoil.cr09.adeg = angle
        self.xfoil.oper()
        return self.xfoil.cr09.cl, self.xfoil.cr09.cd, \
               self.xfoil.cr09.cm, self.xfoil.cl01.lexitflag

    def solveAlphaComplex(self, angle):
        """Compute the flow solution at and angle of attack angle"""

        self.xfoil_cs.cr15.reinf1 = self.re  # Reynolds number
        self.xfoil_cs.cr09.minf1 = self.mach  # Mach Number set
        self.xfoil_cs.ci04.itmax = self.iter  # Iterations Limit Set
        self.xfoil_cs.cr09.adeg = angle
        self.xfoil_cs.oper()
        return self.xfoil_cs.cr09.cl, self.xfoil_cs.cr09.cd, \
               self.xfoil_cs.cr09.cm, self.xfoil_cs.cl01.lexitflag

    def setValue(self, common_block, variable, value):
        """
        Set any xfoil (real) option according to:
        xfoil.<common_block>.<variable> = <value>
        """
        # First check that the attribute the user is trying to set
        # actually exists. These modules work such that you can
        # dynamically create new attributes thus if the varibale
        # doesn't exist it is just created...which is not what we want
        exec("has_attribute = hasattr(self.xfoil.%s,'%s')" % (common_block, variable))
        if has_attribute:  # NOQA: F821
            # Now we can set it
            exec("self.xfoil.%s.%s = value" % (common_block, variable))
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
        exec("has_attribute = hasattr(self.xfoil_cs.%s,'%s')" % (common_block, variable))
        if has_attribute:  # NOQA: F821
            # Now we can set it
            exec("self.xfoil_cs.%s.%s = value" % (common_block, variable))
        else:
            print("There was an error in setValueComplex")
            print("Either %s or %s does not exist" % (common_block, variable))
        # end if
        return

    def xdriver(self):
        cl, cd = self.xfoil.xdriver(self.x, self.y)
        print("Cl=", cl, "Cd=", cd)

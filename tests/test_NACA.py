# Test script to test pyXLIGHT
# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np

# =============================================================================
# Extension modules
# =============================================================================
from pyxlight import pyXLIGHT

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder


class Test_NACA(unittest.TestCase):
    def test_NACA0012(self):
        airfoil = pyXLIGHT.xfoilAnalysis(baseDir + "/" + "naca0012.dat")
        airfoil.re = 100000
        airfoil.mach = 0.0
        airfoil.iter = 100

        print("----------------- Real Results ----------------")
        print("Angle      Cl         Cd         Cm")
        for angle in np.linspace(0, 10, 11):
            cl, cd, cm, lexitflag = airfoil.solveAlpha(angle)
            print("%8f   %8f   %8f   %8f   " % (angle, cl, cd, cm))

        print("----------------- Complex Results ----------------")
        print("Angle      Cl         Cd         Cm")
        for angle in np.linspace(0, 10, 11):
            cl, cd, cm, lexitflag = airfoil.solveAlphaComplex(angle)
            print("%8f   %8f   %8f   %8f   " % (angle, np.real(cl), np.real(cd), np.real(cm)))

    def test_NACA2412(self):
        airfoil = pyXLIGHT.xfoilAnalysis(baseDir + "/" + "naca2412.dat")
        airfoil.re = 100000
        airfoil.mach = 0.0
        airfoil.iter = 100

        print("----------------- Real Results ----------------")
        print("Angle      Cl         Cd         Cm")
        for angle in np.linspace(0, 10, 11):
            cl, cd, cm, lexitflag = airfoil.solveAlpha(angle)
            print("%8f   %8f   %8f   %8f   " % (angle, cl, cd, cm))

        print("----------------- Complex Results ----------------")
        print("Angle      Cl         Cd         Cm")
        for angle in np.linspace(0, 10, 11):
            cl, cd, cm, lexitflag = airfoil.solveAlphaComplex(angle)
            print("%8f   %8f   %8f   %8f   " % (angle, np.real(cl), np.real(cd), np.real(cm)))


if __name__ == "__main__":
    unittest.main()

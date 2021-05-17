# Test script to test pyXLIGHT
# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
from numpy import linspace,real

# =============================================================================
# Extension modules
# =============================================================================
from pyxlight import pyXLIGHT

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder

class Test_NACA(unittest.TestCase):
    def test_NACA0012(self):
        print(os.getcwd())
        airfoil = pyXLIGHT.xfoilAnalysis(baseDir + '/' + 'naca0012.dat')
        airfoil.re = 100000
        airfoil.mach = 0.0
        airfoil.iter = 100

        print('----------------- Real Results ----------------')
        print('Angle      Cl         Cd         Cm')
        for angle in linspace(0,10,11):
            cl,cd,cm,lexitflag = airfoil.solveAlpha(angle)
            print('%8f   %8f   %8f   %8f   '%(angle,cl,cd,cm))

        print('----------------- Complex Results ----------------')
        print('Angle      Cl         Cd         Cm')
        for angle in linspace(0,10,11):
            cl,cd,cm,lexitflag = airfoil.solveAlphaComplex(angle)
            print('%8f   %8f   %8f   %8f   '%(angle,real(cl),real(cd),real(cm)))

    def test_NACA2412(self):
        airfoil = pyXLIGHT.xfoilAnalysis(baseDir + '/' + 'naca2412.dat')
        airfoil.re = 100000
        airfoil.mach = 0.0
        airfoil.iter = 100

        print('----------------- Real Results ----------------')
        print('Angle      Cl         Cd         Cm')
        for angle in linspace(0,10,11):
            cl,cd,cm,lexitflag = airfoil.solveAlpha(angle)
            print('%8f   %8f   %8f   %8f   '%(angle,cl,cd,cm))

        print('----------------- Complex Results ----------------')
        print('Angle      Cl         Cd         Cm')
        for angle in linspace(0,10,11):
            cl,cd,cm,lexitflag = airfoil.solveAlphaComplex(angle)
            print('%8f   %8f   %8f   %8f   '%(angle,real(cl),real(cd),real(cm)))


if __name__ == '__main__':
    unittest.main()
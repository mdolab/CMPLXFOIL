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
from baseclasses import AeroProblem

# =============================================================================
# Extension modules
# =============================================================================
from pyxlight.pyXLIGHT_solver import PYXLIGHT

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder


class Test_NACA(unittest.TestCase):
    def setUp(self):
        self.ap = AeroProblem(
            name="fc", alpha=3, mach=0.2, altitude=1e3, areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cm"]
        )

    def test_NACA0012(self):
        solver = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"))
        funcs = {}
        alphas = np.linspace(-10, 10, 5)
        for i, alpha in enumerate(alphas):
            self.ap.name = f"polar_{alpha}deg"
            self.ap.alpha = alpha
            solver(self.ap)
            solver.evalFunctions(self.ap, funcs)
        true_funcs = {
            "polar_-10.0deg_cd": 0.01101993462628528,
            "polar_-10.0deg_cl": -1.159084142777055,
            "polar_-10.0deg_cm": -0.0011892066512280467,
            "polar_-5.0deg_cd": 0.00662873986434532,
            "polar_-5.0deg_cl": -0.5699935973412055,
            "polar_-5.0deg_cm": -0.0017538538410242423,
            "polar_0.0deg_cd": 0.005119943489125211,
            "polar_0.0deg_cl": -1.5079734336231887e-14,
            "polar_0.0deg_cm": 4.565921255416294e-15,
            "polar_5.0deg_cd": 0.006628892689364835,
            "polar_5.0deg_cl": 0.5700367928924902,
            "polar_5.0deg_cm": 0.0017440066787678468,
            "polar_10.0deg_cd": 0.011018323406864667,
            "polar_10.0deg_cl": 1.1589444225017087,
            "polar_10.0deg_cm": 0.0012188640029523216,
        }
        self.assertDictEqual(funcs, true_funcs)


if __name__ == "__main__":
    unittest.main()

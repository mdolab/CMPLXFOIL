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

externalImportsFailed = False
try:
    from baseclasses import AeroProblem
    from pygeo import DVGeometry, DVGeometryCST
except ImportError:
    externalImportsFailed = True

# =============================================================================
# Extension modules
# =============================================================================
pyxlightImportFailed = False
try:
    from pyxlight import PYXLIGHT
except ImportError:
    pyxlightImportFailed = True

baseDir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder


@unittest.skipIf(externalImportsFailed or pyxlightImportFailed, "Could not import the PYXLIGHT solver")
class TestNACA(unittest.TestCase):
    def setUp(self):
        # Set the random range to use consistent random numbers
        self.rng = np.random.default_rng(7)
        self.ap = AeroProblem(
            name="fc", alpha=3, mach=0.2, altitude=1e3, areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cm"]
        )
        self.CFDSolver = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"))
        self.alphaSequence = np.linspace(-0.5, 0.5, 11)
        self.rng.shuffle(self.alphaSequence)

    def test_cl_solve_default(self):
        """Test that SolveCL works correctly"""
        tol = 1e-5
        for CLTarget in self.alphaSequence:
            # Default settings
            self.CFDSolver.solveCL(self.ap, CLTarget, tol=tol)
            self.assertAlmostEqual(float(self.CFDSolver.xfoil.cr09.cl), CLTarget, delta=tol)

    def test_cl_solve_CLalphaGuess(self):
        """Test that SolveCL works correctly"""
        tol = 1e-5
        for CLTarget in self.alphaSequence:
            # Using secant with CLalphaGuess
            self.CFDSolver.solveCL(self.ap, CLTarget, tol=tol, useNewton=False, CLalphaGuess=0.11)
            self.assertAlmostEqual(float(self.CFDSolver.xfoil.cr09.cl), CLTarget, delta=tol)

    def test_cl_solve_Newton(self):
        """Test that SolveCL works correctly"""
        tol = 1e-5
        for CLTarget in self.alphaSequence:
            # Using Newton
            self.CFDSolver.solveCL(self.ap, CLTarget, tol=tol, useNewton=True)
            self.assertAlmostEqual(float(self.CFDSolver.xfoil.cr09.cl), CLTarget, delta=tol)

    def test_function_values(self):
        funcs = {}
        alphas = np.linspace(-10, 10, 5)
        for _, alpha in enumerate(alphas):
            self.ap.name = f"polar_{alpha}deg"
            self.ap.alpha = alpha
            self.CFDSolver(self.ap)
            self.CFDSolver.evalFunctions(self.ap, funcs)
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
        for key, true_val in true_funcs.items():
            self.assertAlmostEqual(true_val, funcs[key])


@unittest.skipIf(externalImportsFailed or pyxlightImportFailed, "Could not import the PYXLIGHT solver")
class TestTrip(unittest.TestCase):
    def test_trip(self):
        """
        Test that setting the trip changes the result.
        """
        # Set the random range to use consistent random numbers
        self.rng = np.random.default_rng(7)
        evalFuncs = ["cl", "cd", "cm"]
        self.ap = AeroProblem(
            name="fc", alpha=3, mach=0.2, altitude=1e3, areaRef=1.0, chordRef=1.0, evalFuncs=evalFuncs
        )
        self.CFDSolver = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"))
        self.CFDSolverTripped = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"), options={"xTrip": np.array([0.1, 0.1])})
        
        funcs = {}
        self.CFDSolver(self.ap)
        self.CFDSolver.evalFunctions(self.ap, funcs)

        funcsTripped = {}
        self.CFDSolverTripped(self.ap)
        self.CFDSolverTripped.evalFunctions(self.ap, funcsTripped)

        for func in evalFuncs:
            self.assertNotEqual(funcs["fc_" + func], funcsTripped["fc_" + func])


@unittest.skipIf(externalImportsFailed or pyxlightImportFailed, "Could not import the PYXLIGHT solver")
class TestDerivativesFFD(unittest.TestCase):
    def setUp(self):
        # Set the random range to use consistent random numbers
        self.rng = np.random.default_rng(7)

        # Flight conditions
        alpha = 1.5
        mach = 0.3
        Re = 1e6
        T = 288.15  # K

        # Add pyXLIGHT solver
        pyxlightOptions = {
            "writeCoordinates": False,
            "plotAirfoil": False,
            "writeSolution": False,
        }
        self.CFDSolver = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"), options=pyxlightOptions)

        # Add AeroProblem with flight conditions and AoA DV
        self.ap = AeroProblem(
            name="fc",
            alpha=alpha,
            mach=mach,
            reynolds=Re,
            reynoldsLength=1.0,
            T=T,
            areaRef=1.0,
            chordRef=1.0,
            evalFuncs=["cl", "cd", "cm"],
        )
        self.ap.addDV("alpha", value=alpha, lower=0, upper=10.0, scale=1.0)

        # Geometry setup
        FFDFile = os.path.join(baseDir, "naca0012_ffd.xyz")
        self.DVGeo = DVGeometry(FFDFile)
        self.DVGeo.addLocalDV("shape", lower=-0.05, upper=0.05, axis="y", scale=1.0)
        self.CFDSolver.setDVGeo(self.DVGeo)

    def test_alpha_sens(self):
        # Initialize necessary variables
        step = 1e-5  # step size in alpha
        relTol = 5e-3  # acceptable relative error
        alpha = 1.5  # deg
        alphaName = "alpha_" + self.ap.name

        # Evaluate the current functions
        x = {alphaName: alpha}
        self.ap.setDesignVars(x)
        self.CFDSolver(self.ap)
        origFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, origFuncs)

        # Evaluate built in sensitivities with both methods
        funcsSensFD = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensFD, mode="FD")
        self.CFDSolver(self.ap)
        funcsSensCS = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensCS, mode="CS")

        # Perturb alpha
        x[alphaName] += step
        self.ap.setDesignVars(x)
        self.CFDSolver(self.ap)
        pertFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, pertFuncs)

        # Check each function
        for evalFunc in origFuncs.keys():
            checkSensFD = (pertFuncs[evalFunc] - origFuncs[evalFunc]) / step

            # Check evalFunctionsSens's finite difference
            actualSensFD = funcsSensFD[evalFunc][alphaName].item()
            np.testing.assert_allclose(checkSensFD, actualSensFD, rtol=relTol, atol=1e-16)

            # Check evalFunctionsSens's complex step
            actualSensCS = funcsSensCS[evalFunc][alphaName].item()
            np.testing.assert_allclose(checkSensFD, actualSensCS, rtol=relTol, atol=1e-16)

    def test_shape_sens_random(self):
        # Initialize necessary variables
        step = 1e-7  # step size in shape variables
        relTol = 1e-3  # acceptable relative error
        absTol = 1e-3  # acceptable absolute error
        alpha = 1.5  # deg
        alphaName = "alpha_" + self.ap.name
        nShape = 40  # number of shape variables
        shape = self.rng.random(nShape) * 1e-2  # shape variables
        shapeName = "shape"

        # Evaluate the current functions
        x = {alphaName: alpha, shapeName: shape}
        self.ap.setDesignVars(x)
        self.DVGeo.setDesignVars(x)
        self.CFDSolver(self.ap)
        origFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, origFuncs)

        # Evaluate built in sensitivities with both methods
        funcsSensFD = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensFD, mode="FD")
        self.CFDSolver(self.ap)
        funcsSensCS = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensCS, mode="CS")

        # Estimate each partial with finite differences
        checkSensFD = {func: {"shape": np.zeros(nShape)} for func in origFuncs.keys()}  # initialize
        for i in range(nShape):
            x[shapeName][i] += step
            self.DVGeo.setDesignVars(x)
            self.CFDSolver(self.ap)
            pertFuncs = {}
            self.CFDSolver.evalFunctions(self.ap, pertFuncs)
            x[shapeName][i] -= step

            # Compute partial derivative for each function of interest
            for evalFunc in origFuncs.keys():
                checkSensFD[evalFunc][shapeName][i] = (pertFuncs[evalFunc] - origFuncs[evalFunc]) / step

        # Check each function
        for evalFunc in origFuncs.keys():
            # Each partial within each function
            checkVal = checkSensFD[evalFunc][shapeName]
            actualSensFD = funcsSensFD[evalFunc][shapeName]
            actualSensCS = funcsSensCS[evalFunc][shapeName]

            # Check evalFunctionsSens's finite difference
            np.testing.assert_allclose(checkVal, actualSensFD, rtol=relTol, atol=absTol)

            # Check evalFunctionsSens's complex step
            np.testing.assert_allclose(checkVal, actualSensCS, rtol=relTol, atol=absTol)

            # Check evalFunctionsSens's methods against each other
            np.testing.assert_allclose(actualSensFD, actualSensCS, rtol=relTol, atol=absTol)


@unittest.skipIf(externalImportsFailed or pyxlightImportFailed, "Could not import the PYXLIGHT solver")
class TestDerivativesCST(unittest.TestCase):
    def setUp(self):
        # Flight conditions
        self.alpha = 2.
        self.mach = 0.1
        self.Re = 1e7
        self.T = 288.15  # K

        self.nCoeff = 2
        self.DATfile = "n0012BluntTE.dat"

        # Add pyXLIGHT solver
        pyxlightOptions = {
            "writeCoordinates": False,
            "plotAirfoil": False,
            "writeSolution": False,
        }
        self.CFDSolver = PYXLIGHT(os.path.join(baseDir, self.DATfile), options=pyxlightOptions)

        # Add AeroProblem with flight conditions and AoA DV
        self.ap = AeroProblem(
            name="fc",
            alpha=self.alpha,
            mach=self.mach,
            reynolds=self.Re,
            reynoldsLength=1.0,
            T=self.T,
            areaRef=1.0,
            chordRef=1.0,
            evalFuncs=["cl", "cd", "cm"],
        )
        self.ap.addDV("alpha", value=self.alpha, lower=0, upper=10.0, scale=1.0)

        # Geometry setup
        self.CSTName = {"upper": "upper_shape", "lower": "lower_shape"}
        self.CST = {}

        self.CST["upper"] = np.array([0.15, 0.13])
        self.CST["lower"] = np.array([-0.15, -0.1])
        self.upperName = "upper_shape"
        self.lowerName = "lower_shape"

        # This combination of alpha and CST coefficients is likely on a cusp (see multimodality
        # section in the docs), resulting in the FD and CS derivatives mismatching
        # self.alpha = 3.1998091362966057
        # self.CST["upper"] = np.array([0.15335040881179007, 0.1373052558913066])
        # self.CST["lower"] = np.array([-0.15335040881179007, -0.05284956303352384])

        self.DVGeo = DVGeometryCST(os.path.join(baseDir, self.DATfile), debug=False)
        self.DVGeo.addDV(self.CSTName["upper"], dvType="upper", dvNum=self.nCoeff, default=self.CST["upper"])
        self.DVGeo.addDV(self.CSTName["lower"], dvType="lower", dvNum=self.nCoeff, default=self.CST["lower"])
        self.CFDSolver.setDVGeo(self.DVGeo)

    def test_alpha_sens(self):
        # Initialize necessary variables
        step = 1e-6  # step size in alpha
        relTol = 2e-3  # acceptable relative error
        absTol = 1e-5  # acceptable absolute error
        alphaName = "alpha_" + self.ap.name

        # Evaluate the current functions
        x = {alphaName: self.alpha}
        self.ap.setDesignVars(x)
        self.CFDSolver(self.ap)
        origFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, origFuncs)

        # Evaluate built in sensitivities with both methods
        funcsSensFD = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensFD, mode="FD")
        self.CFDSolver(self.ap)
        funcsSensCS = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensCS, mode="CS")

        # Perturb alpha
        x[alphaName] += step
        self.ap.setDesignVars(x)
        self.CFDSolver(self.ap)
        pertFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, pertFuncs)

        # Check each function
        for evalFunc in origFuncs.keys():
            checkSensFD = (pertFuncs[evalFunc] - origFuncs[evalFunc]) / step
            actualSensFD = funcsSensFD[evalFunc][alphaName].item()
            actualSensCS = funcsSensCS[evalFunc][alphaName].item()

            # Check evalFunctionsSens's finite difference
            np.testing.assert_allclose(checkSensFD, actualSensFD, rtol=relTol, atol=absTol)

            # Check evalFunctionsSens's complex step
            np.testing.assert_allclose(checkSensFD, actualSensCS, rtol=relTol, atol=absTol)

    def test_upper_shape_sens(self):
        self._eval_shape_sens("upper")

    def test_lower_shape_sens(self):
        self._eval_shape_sens("lower")

    def _eval_shape_sens(self, surf):
        # Initialize necessary variables
        step = 1e-6  # step size in shape variables
        relTol = 1e-4  # acceptable relative error
        absTol = 1e-4  # acceptable absolute error

        # Evaluate the current functions
        x = {self.CSTName[surf]: self.CST[surf]}
        self.ap.setDesignVars(x)
        self.DVGeo.setDesignVars(x)
        self.CFDSolver(self.ap)
        origFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, origFuncs)

        # Evaluate built in sensitivities with both methods
        funcsSensFD = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensFD, mode="FD")
        self.CFDSolver(self.ap)
        funcsSensCS = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSensCS, mode="CS")

        # Estimate each partial with finite differences
        checkSensFD = {
            func: {self.CSTName[surf]: np.zeros(len(self.CST[surf]))} for func in origFuncs.keys()
        }  # initialize
        for i in range(len(self.CST[surf])):
            x[self.CSTName[surf]][i] += step
            self.DVGeo.setDesignVars(x)
            self.CFDSolver(self.ap)
            pertFuncs = {}
            self.CFDSolver.evalFunctions(self.ap, pertFuncs)
            x[self.CSTName[surf]][i] -= step

            # Compute partial derivative for each function of interest
            for evalFunc in origFuncs.keys():
                checkSensFD[evalFunc][self.CSTName[surf]][i] = (pertFuncs[evalFunc] - origFuncs[evalFunc]) / step

        # Check each function
        for evalFunc in origFuncs.keys():
            # Each partial within each function
            checkVal = checkSensFD[evalFunc][self.CSTName[surf]]
            actualSensFD = funcsSensFD[evalFunc][self.CSTName[surf]]
            actualSensCS = funcsSensCS[evalFunc][self.CSTName[surf]]

            # Check evalFunctionsSens's finite difference
            np.testing.assert_allclose(checkVal, actualSensFD, rtol=relTol, atol=absTol)

            # Check evalFunctionsSens's complex step
            np.testing.assert_allclose(checkVal, actualSensCS, rtol=relTol, atol=absTol)

            # Check evalFunctionsSens's methods against each other
            np.testing.assert_allclose(actualSensFD, actualSensCS, rtol=relTol, atol=absTol)


@unittest.skipIf(externalImportsFailed or pyxlightImportFailed, "Could not import the PYXLIGHT solver")
class TestPlotting(unittest.TestCase):
    """
    Test that the plotting utilities run with no errors.
    """
    def test_plotting(self):
        # Set the random range to use consistent random numbers
        self.rng = np.random.default_rng(13)
        self.ap = AeroProblem(
            name="fc", alpha=3, mach=0.2, altitude=1e3, areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cm"]
        )
        self.CFDSolver = PYXLIGHT(os.path.join(baseDir, "naca0012.dat"))

        # Run an initial case and plot it
        self.CFDSolver(self.ap)
        _, _ = self.CFDSolver.plotAirfoil(showPlot=False)

        # Run another case and update the plot
        self.ap.alpha = 4.
        self.CFDSolver(self.ap)
        _, _ = self.CFDSolver.plotAirfoil(showPlot=False)  # this will call self.CFDSolver.updateAirfoilPlot()

        # Run another case and try calling updateAirfoilPlot directly
        self.ap.alpha = 5.
        self.CFDSolver(self.ap)
        self.CFDSolver.updateAirfoilPlot(pause=False)


if __name__ == "__main__":
    unittest.main()

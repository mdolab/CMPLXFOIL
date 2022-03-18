"""
pyXLIGHT

pyXLIGHT is a wrapper for Mark Drela's Xfoil code. The purpose of this
class is to provide an easy to use wrapper for xfoil for intergration
into other projects. Both real and complex versions of xfoil can be used.

The version uses the MDO Lab's BaseSolver baseclass so it can be used
within the MACH-Aero optimization framework.

Developers:
-----------
- Eytan Adler (EA)
- Alasdair Gray (AG)

History
-------
    v. 1.0 - Initial Class Creation (EA, 2022)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
from baseclasses import BaseSolver
from prefoil.preFoil import readCoordFile

# =============================================================================
# Extension modules
# =============================================================================
from . import MExt
from pyXLIGHT import xfoilAnalysis

class PYXLIGHT(BaseSolver, xfoilAnalysis):
    """
    Solver Class Initialization

    Parameters
    ----------
    fileName : str
        Filename of DAT file to read in
    options : dict, optional
        The user-supplied options dictionary. Options are:
            maxIters
    debug : bool, optional
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging, by default False
    """
    def __init__(self, fileName, options={}, debug=False):
        # Load the compiled module using MExt, allowing multiple imports
        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        time.sleep(0.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil = MExt.MExt("libxlight", curDir, debug=debug)._module
        time.sleep(0.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil_cs = MExt.MExt("libxlight_cs", curDir, debug=debug)._module

        super().__init__("PYXLIGHT", "Panel Method", defaultOptions=self._getDefaultOptions(), options=options)

        # Read in DAT file and create coordinates. DVGeo needs 3D coordinates
        # so keep the third dimension as dummy coordinates.
        self.coords = readCoordFile(fileName)
        self.coords = np.hstack((self.coords, np.zeros((self.coords.shape[0], 1))))
        self.coords0 = self.coords.copy()  # initial coordinates (never changes)

        self.curAP = None

        # Dictionary with dictionary of functions for each aero problem
        self.funcs = {}
        self.functionList = ["cl", "cd", "cm"]  # available functions

    def setDVGeo(self, DVGeo, pointSetKwargs=None):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that PYXLIGHT doe not **strictly** need a
        DVGeometry object, but if optimization with geometric
        changes is desired, then it is required.

        Parameters
        ----------
        DVGeo : A DVGeometry object.
            Object responsible for manipulating the constraints that
            this object is responsible for.
        pointSetKwargs : dict
            Keyword arguments to be passed to the DVGeo addPointSet call.
            Useful for DVGeometryMulti, specifying FFD projection tolerances, etc
        """
        self.DVGeo = DVGeo
        # Get the number of geometry variables
        self.numGeoDV = self.DVGeo.getNDV()

        # Save kwargs for addPointSet
        if pointSetKwargs is not None:
            self.pointSetKwargs = pointSetKwargs

    def setAeroProblem(self, aeroProblem):
        """
        Sets the aeroProblem to by used by PYXLIGHT.

        Parameters
        ----------
        aeroProblem : :class:`AeroProblem <baseclasses.problems.pyAero_problem.AeroProblem>` instance
            Aero problem to set (gives flight conditions)
        """
        ptSetName = "adflow_%s_coords" % aeroProblem.name

        # Tell the user if we are switching aeroProblems
        if self.curAP != aeroProblem:
            self.pp("+" + "-" * 70 + "+")
            self.pp("|  Switching to Aero Problem: %-41s|" % aeroProblem.name)
            self.pp("+" + "-" * 70 + "+")

        # Now check if we have an DVGeo object to deal with:
        if self.DVGeo is not None:
            # DVGeo appeared and we have not embedded points!
            if ptSetName not in self.DVGeo.points:
                self.DVGeo.addPointSet(self.coords0, ptSetName, **self.pointSetKwargs)

            # Check if our point-set is up to date:
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                coords = self.DVGeo.update(ptSetName, config=aeroProblem.name)

                self.setCoordinates(coords)

        self.curAP = aeroProblem

    def setCoordinates(self, coords):
        """
        Update the airfoil coordinates and associated point sets.

        Parameters
        ----------
        coords : ndarray
            New airfoil coordinates (either 2 or 3 columns)
        """
        self.coords[:, :2] = coords[:, :2].copy()
        super().setCoordinates(self.coords[:, 0], self.coords[:, 1])

    def getCoordinates(self):
        """
        Return the current airfoil coordinates

        Returns
        -------
        coords : ndarray
            Airfoil coordinates with each column being (x, y, z)
            where z is a dummy value.
        """
        return self.coords.copy()

    def __call__(self, aeroProblem):
        """
        Evaluate XFOIL with the current coordinates and flight conditions (from aeroProblem).

        Parameters
        ----------
        aeroProblem : :class:`AeroProblem <baseclasses.problems.pyAero_problem.AeroProblem>` instance
            Aero problem to set (gives flight conditions)
        """
        self.setAeroProblem(aeroProblem)

        # Set flight condition and options
        self.xfoil.cr15.reinf1 = aeroProblem.reynolds  # Reynolds number
        self.xfoil.cr09.minf1 = aeroProblem.mach  # Mach Number set
        self.xfoil.cr09.adeg = aeroProblem.alpha
        self.xfoil.ci04.itmax = self.getOption("maxIters")  # Iterations Limit Set

        # Call XFOIL
        self.xfoil.oper()

        # Store results in dictionary for current aero problem
        self.funcs[aeroProblem.name] = {
            "cl": self.xfoil.cr09.cl,
            "cd": self.xfoil.cr09.cd,
            "cm": self.xfoil.cr09.cm,
        }

        # CHeck for failure
        self.curAP.solveFailed = self.curAP.fatalFail = self.xfoil.cl01.lexitflag[0] != 0

    def checkSolutionFailure(self, aeroProblem, funcs):
        """Take in a an aeroProblem and check for failure.

        Then append the fail flag in funcs. Information regarding whether or not the last analysis with the aeroProblem
        was sucessful is included. This information is included as "funcs['fail']". If the 'fail' entry already exits in
        the dictionary the following operation is performed:

        funcs['fail'] = funcs['fail'] or <did this problem fail>

        In other words, if any one problem fails, the funcs['fail'] entry will be True. This information can then be
        used directly in multiPointSparse. For direct interface with pyOptSparse the fail flag needs to be returned
        separately from the funcs.

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to to get the solution for
        funcs : dict
            Dictionary into which the functions are saved.
        """

        self.setAeroProblem(aeroProblem)
        # We also add the fail flag into the funcs dictionary. If fail is already there, we just logically 'or' what was
        # there. Otherwise we add a new entry.
        failFlag = self.curAP.solveFailed or self.curAP.fatalFail
        if "fail" in funcs:
            funcs["fail"] = funcs["fail"] or failFlag
        else:
            funcs["fail"] = failFlag

    def evalFunctions(self, aeroProblem, funcs, evalFuncs=None, ignoreMissing=False):
        """
        This is the main routine for returning useful information from PYXLIGHT.
        The functions corresponding to the strings in ``evalFuncs`` are evaluated
        and updated into the provided dictionary.

        Parameters
        ----------
        aeroProblem : :class:`AeroProblem <baseclasses.problems.pyAero_problem.AeroProblem>` instance
            Aero problem from which to pull evalFuncs and flight conditions.
        funcs : dict
            Dictionary into which the functions are saved.
        evalFuncs : iterable object containing strings
            If not none, use these functions to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid function. Please use
            this option with caution.
        """
        self.setAeroProblem(aeroProblem)
        if evalFuncs is None:
            evalFuncs = sorted(list(self.curAP.evalFuncs))
        else:o
            evalFuncs = sorted(list(evalFuncs))

        # Make the functions lower case
        evalFuncs = [s.lower() for s in evalFuncs]

        returnFuncs = {}
        for f in evalFuncs:
            # If it's not in the list of available functions
            if f not in self.functionList:
                # Either throw an error (if requested) or skip it
                if not ignoreMissing:
                    raise ValueError(f"Supplied function \"{f}\" is not in the available functions {self.functionList}.")
            else:
                returnFuncs[aeroProblem.name + "_" + f] = self.funcs[aeroProblem.name][f]
        funcs.update(returnFuncs)

    @staticmethod
    def _getDefaultOptions():
        return {
            "maxIters": [int, 100],  # maximum iterations for XFOIL solver
        }

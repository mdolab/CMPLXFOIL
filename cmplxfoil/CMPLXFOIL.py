"""
CMPLXFOIL

CMPLXFOIL is a wrapper for Mark Drela's XFOIL code. The purpose of this
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
    v. 1.0 - Initial Class Creation (EA + AG, 2022)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
from copy import copy, deepcopy
import pickle as pkl
from collections.abc import Iterable
import warnings

# =============================================================================
# Other Python modules
# =============================================================================
import numpy as np
from baseclasses import BaseSolver
from . import MExt

# Handle plotting imports
plt = None
try:
    from matplotlib import pyplot
    import matplotlib.lines as mpl_lines

    plt = pyplot
except ImportError:
    pass

# If the matplotlib import worked, try niceplots
if plt:
    try:
        import niceplots

        plt.style.use(niceplots.get_style())
        colors = niceplots.get_colors()
        color = colors["Blue"]
        cpUpColor = colors["Blue"]
        cpLowColor = colors["Red"]
    except ImportError:
        color = "b"
        cpUpColor = "b"
        cpLowColor = "r"


class CMPLXFOIL(BaseSolver):
    """
    CMPLXFOIL Class Initialization

    Parameters
    ----------
    fileName : str
        Filename of DAT file to read in
    options : dict of option-value pairs, optional
        Options for the solver. Available options can be found
        in the Options section of the documentation or the
        options.yaml file in the docs directory.
    debug : bool, optional
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging, by default False
    """

    def __init__(self, fileName, options={}, debug=False):
        # Load the compiled module using MExt, allowing multiple imports
        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        time.sleep(0.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil = MExt.MExt("libcmplxfoil", curDir, debug=debug)._module
        time.sleep(0.1)  # BOTH of these sleeps are necessary for some reason!
        self.xfoil_cs = MExt.MExt("libcmplxfoil_cs", curDir, debug=debug)._module

        super().__init__("CMPLXFOIL", "Panel Method", defaultOptions=self._getDefaultOptions(), options=options)

        if self.getOption("xTrip").shape != (2,):
            raise ValueError("xTrip option shape must be (2,)")

        # Read in DAT file and create coordinates. DVGeo needs 3D coordinates
        # so keep the third dimension as dummy coordinates.
        self.DATFileName = fileName
        self.coords = self._readDat(fileName)
        self.coords = np.hstack((self.coords, np.zeros((self.coords.shape[0], 1))))
        self.coords0 = self.coords.copy()  # initial coordinates (never changes)
        self.setCoordinates(self.coords0)  # set the initial coordinates
        self.setCoordinatesComplex(self.coords0)  # set the initial complex coordinates
        self.pointSetKwargs = {}

        self.curAP = None
        self.DVGeo = None

        # Dictionary with dictionary of functions for each aero problem
        self.funcs = {}
        self.funcsComplex = {}
        self.functionList = ["cl", "cd", "cm", "kscpmin"]  # available functions

        # When the XFOIL solver is called, slice data is saved (key is the current
        # AeroProblem name). In the associated value is a dictionary containing
        #       "cp_visc_upper": viscous CP on the airfoil's upper surface
        #       "cp_invisc_upper": inviscid CP on the airfoil's upper surface
        #       "x_upper": x coordinates of the upper surface CP data
        #       "y_upper": y coordinates of the upper surface CP data
        #       "cp_visc_lower": viscous CP on the airfoil's lower surface
        #       "cp_invisc_lower": inviscid CP on the airfoil's lower surface
        #       "x_lower": x coordinates of the lower surface CP data
        #       "y_lower": y coordinates of the lower surface CP data
        #       "cf_upper": skin friction coefficient on the upper surface
        #       "x_cf_upper": x coordinates of upper surface skin friction coefficient
        #       "y_cf_upper": y coordinates of upper surface skin friction coefficient
        #       "cf_lower": skin friction coefficient on the lower surface
        #       "x_cf_lower": x coordinates of lower surface skin friction coefficient
        #       "y_cf_lower": y coordinates of lower surface skin friction coefficient
        self.sliceData = {}
        self.sliceDataComplex = {}

        # Possible AeroProblem design variables (only alpha for CMPLXFOIL)
        self.possibleAeroDVs = ["alpha"]

        # Figure and axes used by self.plotAirfoil and self._updateAirfoilPlot
        self.airfoilFig = None
        self.airfoilAxs = None
        self.CPlim = None  # y limits on the CP plot
        self.CFlim = None  # y limits on the CF plot

    def setDVGeo(self, DVGeo, pointSetKwargs=None):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that CMPLXFOIL does not **strictly** need a
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
        Sets the aeroProblem to by used by CMPLXFOIL.

        Parameters
        ----------
        aeroProblem : :class:`AeroProblem <baseclasses.problems.pyAero_problem.AeroProblem>` instance
            Aero problem to set (gives flight conditions)
        """
        ptSetName = "cmplxfoil_%s_coords" % aeroProblem.name

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
                aeroProblem.ptSetName = ptSetName

            # Check if our point-set is up to date:
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                coords = self.DVGeo.update(ptSetName, config=aeroProblem.name)

                self.setCoordinates(coords)

        # Initialize the callCounter if it's not already an attribute
        try:
            _ = aeroProblem.callCounter
        except AttributeError:
            aeroProblem.callCounter = 0

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
        self._setCoordinates(self.coords[:, 0], self.coords[:, 1])

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

    def setCoordinatesComplex(self, coords):
        """
        Update the complex airfoil coordinates and associated point sets.

        Parameters
        ----------
        coords : ndarray
            New airfoil coordinates (either 2 or 3 columns)
        """
        self._setCoordinates(coords[:, 0].astype(complex), coords[:, 1].astype(complex), setComplex=True)

    def __call__(self, aeroProblem, useComplex=False, deriv=False):
        """
        Evaluate XFOIL with the current coordinates and flight conditions (from aeroProblem).

        Parameters
        ----------
        aeroProblem : :class:`AeroProblem <baseclasses.problems.pyAero_problem.AeroProblem>` instance
            Aero problem to set (gives flight conditions)
        useComplex : bool
            Run XFOIL in complex mode
        deriv : bool
            Set to True if the call is associated with the derivative calculation. callCounter will be incremented and
            writeSolution called only if this is False
        """
        if useComplex:
            dtype = complex
            xfoil = self.xfoil_cs
            funcs = self.funcsComplex
            sliceData = self.sliceDataComplex
            xfoil.cl01.printconv = self.getOption("printComplexConvergence")
        else:
            dtype = float
            xfoil = self.xfoil
            funcs = self.funcs
            sliceData = self.sliceData
            xfoil.cl01.printconv = self.getOption("printRealConvergence")

        self.setAeroProblem(aeroProblem)

        # Set flight condition and options
        xfoil.cr15.reinf1 = aeroProblem.re  # Reynolds number per unit length
        xfoil.cr09.minf1 = aeroProblem.mach  # Mach number
        xfoil.cr09.adeg = aeroProblem.alpha  # Angle of attack
        xfoil.ci04.itmax = self.getOption("maxIters")  # Iterations limit
        if not np.any(np.isnan(self.getOption("xTrip"))):  # nan is default to not set, otherwise set it
            xfoil.cr15.xstrip = self.getOption("xTrip")

        # Set nCrit (The Fortran variable is acrit)
        # lvconv needs to be set to False to update internal variables
        xfoil.cr15.acrit = self.getOption("nCrit")
        xfoil.cl01.lvconv = False

        # Call XFOIL
        xfoil.oper()

        # Store results in dictionary for current aero problem
        funcs[aeroProblem.name] = {
            "cl": dtype(xfoil.cr09.cl),
            "cd": dtype(xfoil.cr09.cd),
            "cm": dtype(xfoil.cr09.cm),
        }

        # Pull out and process the pressure and skin friction coefficient data
        cpv = xfoil.cr04.cpv  # viscous
        cpi = xfoil.cr04.cpi  # inviscid
        x = xfoil.cr05.x
        y = xfoil.cr05.y
        end_foil_idx = np.argmax(x > self.coords[0, 0]) + 1  # XFOIL includes wake panels, which we don't want
        idx_lower_start = end_foil_idx // 2  # first half of data is upper surface

        idxUpper = np.argwhere(xfoil.cr15.tau[:, 0] != 0).flatten()
        idxLower = np.argwhere(xfoil.cr15.tau[:, 1] != 0).flatten()
        iPanUpper = xfoil.ci05.ipan[idxUpper, 0] - 1  # FORTRAN uses 1-based indexing, need to adjust
        iPanLower = xfoil.ci05.ipan[idxLower, 1] - 1  # FORTRAN uses 1-based indexing, need to adjust
        tauUpper = xfoil.cr15.tau[idxUpper, 0].copy().astype(dtype)
        tauLower = xfoil.cr15.tau[idxLower, 1].copy().astype(dtype)
        uInf = xfoil.cr09.qinf
        xCfUpper = x[iPanUpper].copy().astype(dtype)
        yCfUpper = y[iPanUpper].copy().astype(dtype)
        xCfLower = x[iPanLower].copy().astype(dtype)
        yCfLower = y[iPanLower].copy().astype(dtype)

        sliceData[aeroProblem.name] = {
            "cp_visc_upper": cpv[:idx_lower_start].copy().astype(dtype),
            "cp_invisc_upper": cpi[:idx_lower_start].copy().astype(dtype),
            "x_upper": x[:idx_lower_start].copy().astype(dtype),
            "y_upper": y[:idx_lower_start].copy().astype(dtype),
            "cp_visc_lower": cpv[idx_lower_start:end_foil_idx].copy().astype(dtype),
            "cp_invisc_lower": cpi[idx_lower_start:end_foil_idx].copy().astype(dtype),
            "x_lower": x[idx_lower_start:end_foil_idx].copy().astype(dtype),
            "y_lower": y[idx_lower_start:end_foil_idx].copy().astype(dtype),
            "cf_upper": tauUpper / (0.5 * uInf**2),
            "x_cf_upper": xCfUpper,
            "y_cf_upper": yCfUpper,
            "cf_lower": tauLower / (0.5 * uInf**2),
            "x_cf_lower": xCfLower,
            "y_cf_lower": yCfLower,
        }

        # Check if kscpmin is requested, if so, then compute it
        if "kscpmin" in self.curAP.evalFuncs:
            cpAll = np.concatenate((sliceData[aeroProblem.name]["cp_visc_upper"], sliceData[aeroProblem.name]["cp_visc_lower"]))
            kscpmin = -self.computeKSMax(-cpAll, rho=self.getOption("rhoKS"), printOK=False)

            funcs[aeroProblem.name]["kscpmin"] = dtype(kscpmin)

        # Check for failure
        self.curAP.solveFailed = self.curAP.fatalFail = xfoil.cl01.lexitflag != 0 or xfoil.cl01.lvconv == 0

        # If not a derivative call, increment callCounter
        if not deriv:
            self.curAP.callCounter += 1

        # Write solution files if desired
        if not deriv and self.getOption("writeSolution"):
            self.writeSolution()

    def computeKSMax(self, g, rho=10000.0, printOK=True):
        """
        Compute a smooth approximation to the maximum of a set of values
        using Kreisselmeier--Steinhauser aggregation.

        Parameters
        ----------
        g : 1d array
            Values of which to approximate the maximum
        rho : float, optional
            KS Weight parameter, larger values give a closer but less smooth
            approximation of the maximum, by default 10000.0

        Returns
        -------
        ksmax : float
            The KS aggregated value
        """

        maxg = np.max(g)
        ksmax = maxg + 1.0 / rho * np.log(np.sum(np.exp(rho * (g - maxg))))

        if printOK:
            print(f"true max: {maxg} \nks max:   {ksmax}")

        return ksmax

    def solveCL(
        self,
        aeroProblem,
        CLStar,
        alpha0=None,
        alphaBound=None,
        delta=0.5,
        tol=1e-3,
        CLalphaGuess=None,
        maxIter=20,
        useNewton=False,
    ):
        """Find the angle of attack that gives a target lift coefficient.

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to solve
        CLStar : float
            The desired CL
        alpha0 : float, optional
            Initial guess for secant search (deg). If None, use the value in the aeroProblem, by default None
        alphaBound : float, tuple, list, optional
            Bounds for angle of attack, if scalar then value is treated as a +- bound, by default None, in which case
            limit is +/-15 deg
        delta : float, optional
            Initial step direction for secant search, by default 0.5
        tol : float, optional
            Desired tolerance for CL, by default 1e-3
        CLalphaGuess : float, optional
            The user can provide an estimate for the lift curve slope in order to accelerate convergence. If the user
            supplies a value to this option, it will not use the delta value anymore to select the angle of attack of
            the second run. The value should be in 1/deg., by default None
        maxIter : int, optional
            Maximum number of iterations, by default 20
        useNewton : bool, optional
            If True, Newton's method will be used where the dCL/dAlpha is computed using complex-step, otherwise the
            secant method is used, by default False

        Returns
        -------
        None, but the correct alpha is stored in the aeroProblem
        """
        self.setAeroProblem(aeroProblem)
        if alpha0 is not None:
            aeroProblem.alpha = alpha0
        else:
            alpha0 = aeroProblem.alpha

        if alphaBound is None:
            alphaBound = (-15, 15)
        elif isinstance(alphaBound, (int, float)):
            alphaBound = (-alphaBound, alphaBound)
        elif isinstance(alphaBound, Iterable):
            alphaBound = (alphaBound[0], alphaBound[1])
        else:
            raise ValueError(
                f'Supplied alphaBound value "{alphaBound}" is not the correct type, must be a scalar or array-like value.'
            )

        dCLdAlpha = CLalphaGuess
        resPrev = None
        alphaPrev = None
        hasConverged = False
        maxRes = 1e2

        for _ in range(maxIter):
            # Call the solver and compute the "residual"
            self.__call__(aeroProblem, deriv=True)
            failCheck = {}
            self.checkSolutionFailure(aeroProblem, failCheck)
            cl = float(self.funcs[aeroProblem.name]["cl"])
            print(f"Alpha: {aeroProblem.alpha:>6.3f}, CL: {cl:>7.6f}")
            res = cl - CLStar

            # Check for convergence and divergence
            hasConverged = abs(res) < tol and not failCheck["fail"]
            if hasConverged:
                print("Converged!")
                break
            elif abs(res > maxRes):
                warnings.warn(f"solveCL diverged (CL = {cl:.2f})", stacklevel=2)
                break

            # If not converged or diverged, compute the next alpha
            if not failCheck["fail"]:
                if useNewton:
                    aeroProblem.alpha += 1e-200 * 1j
                    self.__call__(aeroProblem, useComplex=True, deriv=True)
                    dCLdAlpha = np.imag(self.funcsComplex[aeroProblem.name]["cl"]) * 1e200
                    aeroProblem.alpha = np.real(aeroProblem.alpha)
                else:
                    # If using secant, we can only compute dCLdAlpha from the second iteration onwards
                    if resPrev is not None:
                        dCLdAlpha = (res - resPrev) / (aeroProblem.alpha - alphaPrev)
                resPrev = res
                alphaPrev = aeroProblem.alpha

            # Update the alpha either using dCLdAlpha or delta, or backtracking if something went wrong
            if dCLdAlpha == 0.0 or failCheck["fail"]:
                aeroProblem.alpha *= 0.9
            elif dCLdAlpha is not None:
                aeroProblem.alpha = np.clip(aeroProblem.alpha - res / dCLdAlpha, alphaBound[0], alphaBound[1])
            else:
                aeroProblem.alpha = aeroProblem.alpha + delta
        if not hasConverged:
            print("Did not converge :(")
        return hasConverged

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """This is a generic shell function that potentially writes the various output files. The intent is that the
        user or calling program can call this file and CMPLXFOIL write all the files that the user has defined. It is
        recommended that this function is used along with the associated logical flags in the options to determine the
        desired writing procedure.

        Parameters
        ----------
        outputDir : str
            Use the supplied output directory
        baseName : str
            Use this supplied string for the base filename. Typically only used from an external solver.
        number : int
            Use the user supplied number to index solution. Again, only typically used from an external solver.
        """
        if outputDir is None:
            outputDir = self.getOption("outputDirectory")

        if baseName is None:
            baseName = self.curAP.name

        # Add a number to the filename, either from the user or from the current callCounter
        if number is not None:
            baseName = baseName + "_%3.3d" % number
        else:
            if self.getOption("numberSolutions"):
                baseName = baseName + "_%3.3d" % self.curAP.callCounter

        # Join to get the actual filename root
        base = os.path.join(outputDir, baseName)

        if self.getOption("writeCoordinates"):
            self.writeCoordinates(base)

        if self.getOption("writeSliceFile"):
            self.writeSlice(base)

        if self.getOption("plotAirfoil"):
            self.plotAirfoil()

    def writeCoordinates(self, fileName):
        """Write dat file with the current coordinates.

        Parameters
        ----------
        fileName : str
            File name for saved dat file (".dat" will be automatically appended).
        """
        fileName += ".dat"
        with open(fileName, "w") as f:
            for i in range(self.coords.shape[0]):
                f.write(str(round(self.coords[i, 0], 12)) + "\t\t" + str(round(self.coords[i, 1], 12)) + "\n")

    def writeSlice(self, fileName):
        """Write pickle file containing the sliceData dictionary. The data can be
        accessed using the AeroProblem name as the key. Within that is a dictionary containing

        * Pressure coefficient data on the upper surface

            * ``"cp_visc_upper"``: viscous CP on the airfoil's upper surface
            * ``"cp_invisc_upper"``: inviscid CP on the airfoil's upper surface
            * ``"x_upper"``: x coordinates of the upper surface CP data
            * ``"y_upper"``: y coordinates of the upper surface CP data

        * Pressure coefficient data on the lower surface

            * ``"cp_visc_lower"``: viscous CP on the airfoil's lower surface
            * ``"cp_invisc_lower"``: inviscid CP on the airfoil's lower surface
            * ``"x_lower"``: x coordinates of the lower surface CP data
            * ``"y_lower"``: y coordinates of the lower surface CP data

        * Skin friction coefficient data on the upper surface

            * ``"cf_upper"``: skin friction coefficient on the upper surface
            * ``"x_cf_upper"``: x coordinates of upper surface skin friction coefficient
            * ``"y_cf_upper"``: y coordinates of upper surface skin friction coefficient

        * Skin friction coefficient data on the lower surface

            * ``"cf_lower"``: skin friction coefficient on the lower surface
            * ``"x_cf_lower"``: x coordinates of lower surface skin friction coefficient
            * ``"y_cf_lower"``: y coordinates of lower surface skin friction coefficient

        Parameters
        ----------
        fileName : str
            File name for saved pkl file (".pkl" will be automatically appended).
        """
        fileName += ".pkl"
        with open(fileName, "wb") as f:
            pkl.dump(self.sliceData, f, protocol=pkl.HIGHEST_PROTOCOL)

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
            The aerodynamic problem to get the solution for
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

    def checkAdjointFailure(self, aeroProblem, funcsSens):
        """
        Pass through to checkSolutionFailure to maintain the same interface as ADflow.

        This checks if the primal solve fails and can be called when the sensitivity
        is being evaluated (either through FD or CS).

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to to get the solution for
        funcsSens : dict
            Dictionary into which the functions are saved.
        """
        self.checkSolutionFailure(aeroProblem, funcsSens)

    def evalFunctions(self, aeroProblem, funcs, evalFuncs=None, ignoreMissing=False):
        """
        This is the main routine for returning useful information from CMPLXFOIL.
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
        else:
            evalFuncs = sorted(list(evalFuncs))

        # Make the functions lower case
        evalFuncs = [s.lower() for s in evalFuncs]

        returnFuncs = {}
        for f in evalFuncs:
            # If it's not in the list of available functions
            if f not in self.functionList:
                # Either throw an error (if requested) or skip it
                if not ignoreMissing:
                    raise ValueError(f'Supplied function "{f}" is not in the available functions {self.functionList}.')
            else:
                returnFuncs[aeroProblem.name + "_" + f] = self.funcs[aeroProblem.name][f]
        funcs.update(returnFuncs)

        # Set the AeroProblem's function names
        for name in evalFuncs:
            self.curAP.funcNames[name] = self.curAP.name + "_" + name

    def evalFunctionsSens(self, aeroProblem, funcsSens, evalFuncs=None, mode="CS", h=None):
        """
        Evaluate the sensitivity of the desired functions given in
        iterable object, 'evalFuncs' and add them to the dictionary
        'funcSens'. The keys in the 'funcsSens' dictionary will be have an
        ``<ap.name>_`` prepended to them.

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the function derivatives are saved
        evalFuncs : iterable object containing strings
            The additional functions the user wants returned that are
            not already defined in the aeroProblem
        mode : str ["FD" or "CS"]
            Specifies how the jacobian vector products will be computed
        h : float
            Step sized used when the mode is "FD" or "CS" (must be complex
            if mode = "CS")
        """
        # This is the one and only gateway to the getting derivatives
        # out of xfoil. If you want a derivative, it should come from
        # here. Thank you.
        self.setAeroProblem(aeroProblem)
        if evalFuncs is None:
            evalFuncs = sorted(list(self.curAP.evalFuncs))
        else:
            evalFuncs = sorted(list(evalFuncs))

        # Make the functions lower case
        evalFuncs = [s.lower() for s in evalFuncs]

        # Get design variables
        DVs = self.DVGeo.getValues()
        for dv in self.curAP.DVs.values():
            DVs[dv.key] = np.atleast_1d(dv.value)

        # Preallocate the funcsSens dictionary with zeros for the desired sensitivities
        for f in evalFuncs:
            funcsSens[self.curAP.name + "_" + f] = {}
            for dvName, dvVal in DVs.items():
                if isinstance(dvVal, np.ndarray):
                    funcsSens[self.curAP.name + "_" + f][dvName] = np.zeros(dvVal.shape, dtype=float)
                else:
                    funcsSens[self.curAP.name + "_" + f][dvName] = np.zeros(1)

        # Loop over design variables and compute derivatives of each
        for dvName, dvVal in DVs.items():
            lenDV = 1 if not isinstance(dvVal, np.ndarray) else len(dvVal)
            for i in range(lenDV):
                # Compute the seed for the finite difference/complex step
                seed = {dvName: np.zeros(dvVal.shape, dtype=float)}
                seed[dvName][i] = 1.0

                # Compute the sensitivity
                sens = self.computeJacobianVectorProductFwd(xDvDot=seed, mode=mode, h=h)
                for f in evalFuncs:
                    funcsSens[self.curAP.name + "_" + f][dvName][i] = sens[f]

                # Check that the solution converged
                self.checkSolutionFailure(self.curAP, funcsSens)

        # Append "_" + (current aero problem name) to any design variables associated with the aero problem
        for f in evalFuncs:
            func = self.curAP.name + "_" + f
            for dvName in DVs.keys():
                if dvName in self.possibleAeroDVs and dvName in funcsSens[func]:
                    funcsSens[func][dvName + "_" + self.curAP.name] = funcsSens[func][dvName]
                    del funcsSens[func][dvName]

    def computeJacobianVectorProductFwd(self, xDvDot=None, xSDot=None, mode="CS", h=None):
        """This the main Python gateway for producing forward mode jacobian
        vector products. They are computed using either the complex step or finite
        difference method. This function is not generally called by the user but
        rather internally or from another solver. A DVGeo object must be set
        for this routine.

        Parameters
        ----------
        xDvDot : dict
            Perturbation on the design variables
        xSDot : numpy array
            Perturbation on the surface
        mode : str ["FD" or "CS"]
            Specifies how the jacobian vector products will be computed
        h : float
            Step sized used when the mode is "FD" or "CS" (must be complex
            if mode = "CS"), by default 1e-6 for FD and 1e-200j for CS

        Returns
        -------
        dict
            Jacobian vector product of evalFuncs given perturbation
        """
        if xDvDot is None and xSDot is None:
            raise ValueError("xDvDot and xSDot cannot both be None")

        if mode not in ["FD", "CS"]:
            raise ValueError(f'Jacobian vector product mode "{mode}" invalid. Must be either "FD" or "CS"')

        if self.DVGeo is None:
            raise ValueError("DVGeo object must be added with setDVGeo before calling computeJacobianVectorProductFwd")

        geoDVs = list(self.DVGeo.getValues().keys())
        possibleDVs = self.possibleAeroDVs + geoDVs
        for DV in xDvDot.keys():
            if DV not in possibleDVs:
                raise ValueError(f'Perturbed design variable "{DV}" is not valid')

        if h is None:
            if mode == "FD":
                h = 1e-6
            elif mode == "CS":
                h = 1e-200j

        if mode == "FD":
            orig_funcs = deepcopy(self.funcs[self.curAP.name])

        # Process the Xs perturbation
        if xSDot is None:
            xsdot = np.zeros_like(self.coords0)
        else:
            xsdot = xSDot

        # For the geometric xDvDot perturbation we accumulate into the
        # already existing (and possibly nonzero) xsdot and xvdot
        geoxDvDot = {k: xDvDot[k] for k in geoDVs if k in xDvDot}
        if xDvDot is not None and self.DVGeo is not None:
            xsdot += self.DVGeo.totalSensitivityProd(geoxDvDot, self.curAP.ptSetName, config=self.curAP.name).reshape(
                xsdot.shape
            )

        # Perturb the coordinates
        orig_coords = self.getCoordinates()
        if mode == "FD":
            self.setCoordinates(orig_coords + xsdot * h)
        else:
            self.setCoordinatesComplex(orig_coords + xsdot * h)

        orig_alpha = copy(self.curAP.alpha)
        if "alpha" in xDvDot.keys():
            self.curAP.alpha += xDvDot["alpha"] * h

        self.__call__(self.curAP, useComplex=mode == "CS", deriv=True)

        # Compute the Jacobian vector products
        jacVecProd = {}
        for f in self.functionList:
            if mode == "FD":
                jacVecProd[f] = (self.funcs[self.curAP.name][f] - orig_funcs[f]) / h
            else:
                jacVecProd[f] = np.imag(self.funcsComplex[self.curAP.name][f]) / np.imag(h)

        # Reset the perturbed variables
        self.curAP.alpha = orig_alpha
        if mode == "FD":
            self.setCoordinates(orig_coords)
            self.funcs[self.curAP.name] = orig_funcs
        else:
            self.setCoordinatesComplex(orig_coords)

        return jacVecProd

    def getTriangulatedMeshSurface(self, offsetDist=1.0):
        """
        This function returns a pyGeo surface. The intent is
        to use this for DVConstraints.

        .. note::
            This method requires the pyGeo library

        Parameters
        ----------
        offsetDist : float
            Distance to extrude airfoil (same units as airfoil coordinates)

        Returns
        -------
        pyGeo surface
            Extruded airfoil surface
        """
        try:
            from pygeo.pyGeo import pyGeo
        except ImportError as e:
            raise ImportError("pygeo is required to use getTriangulatedMeshSurface") from e

        airfoil_list = [self.DATFileName] * 2
        naf = len(airfoil_list)

        # Airfoil leading edge positions
        x = [0.0, 0.0]
        y = [0.0, 0.0]
        z = [0.0, offsetDist]
        offset = np.zeros((naf, 2))  # x-y offset applied to airfoil position before scaling

        # Airfoil rotations
        rot_x = [0.0, 0.0]
        rot_y = [0.0, 0.0]
        rot_z = [0.0, 0.0]

        # Airfoil scaling
        scale = [1.0, 1.0]  # scaling factor on chord lengths

        # Find the trailing edge thickness
        teThickness = self.coords[0, 1] - self.coords[-1, 1]

        return pyGeo(
            "liftingSurface",
            xsections=airfoil_list,
            scale=scale,
            offset=offset,
            x=x,
            y=y,
            z=z,
            rotX=rot_x,
            rotY=rot_y,
            rotZ=rot_z,
            bluntTe=True,
            squareTeTip=True,
            teHeight=teThickness,
        )

    def plotAirfoil(self, fileName=None, showPlot=True):
        """
        Plots the current airfoil and returns the figure.

        Parameters
        ----------
        fileName : str, optional
            FileName to save to, if none specified it will
            show the plot with plt.show()
        showPlot : bool, optional
            Pop open the plot, by default True

        Returns
        -------
        matplotlib figure
            Figure with airfoil plotted to it
        list of matplotlib axes
            List of matplotlib axes for CP and airfoil plots (in that order)
        """
        if plt is None:
            raise ImportError("matplotlib is required to use the plotting tools")

        if self.airfoilFig is None:
            # Get data to plot
            x = self.coords[:, 0]
            y = self.coords[:, 1]
            CPUpper = self.sliceData[self.curAP.name]["cp_visc_upper"]
            CPUpperInvisc = self.sliceData[self.curAP.name]["cp_invisc_upper"]
            xUpper = self.sliceData[self.curAP.name]["x_upper"]
            CPLower = self.sliceData[self.curAP.name]["cp_visc_lower"]
            CPLowerInvisc = self.sliceData[self.curAP.name]["cp_invisc_lower"]
            xLower = self.sliceData[self.curAP.name]["x_lower"]
            CFUpper = self.sliceData[self.curAP.name]["cf_upper"]
            xCFUpper = self.sliceData[self.curAP.name]["x_cf_upper"]
            CFLower = self.sliceData[self.curAP.name]["cf_lower"]
            xCFLower = self.sliceData[self.curAP.name]["x_cf_lower"]

            # Inverted CP axis
            stackedCP = np.hstack((CPUpper, CPUpperInvisc, CPLower, CPLowerInvisc))
            stackedCP = stackedCP[np.isfinite(stackedCP)]
            stackedCF = np.hstack((CFUpper, CFLower))
            stackedCF = stackedCF[np.isfinite(stackedCF)]
            self.CPlim = [max(stackedCP) + 0.2, min(stackedCP) - 0.2]
            self.CFlim = [min(stackedCF) - 0.002, max(stackedCF) + 0.002]
            self.xlimFoil = [min(x) - 0.01, max(x) + 0.01]
            self.ylimFoil = [min(y) - 0.01, max(y) + 0.01]

            fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=[10, 13])
            iAxsCP = 0
            iAxsCF = 1
            iAxsFoil = 2
            plt.ion()
            if showPlot:
                plt.show()

            # Plot the CP on the upper axis
            axs[iAxsCP].plot(xUpper, CPUpper, color="k", zorder=-1, alpha=0.15)
            axs[iAxsCP].plot(xLower, CPLower, color="k", zorder=-1, alpha=0.15)
            axs[iAxsCP].plot(xUpper, CPUpper, color=cpUpColor)
            axs[iAxsCP].plot(xLower, CPLower, color=cpLowColor)
            axs[iAxsCP].plot(xUpper, CPUpperInvisc, "--", color=cpUpColor, linewidth=1.0)
            axs[iAxsCP].plot(xLower, CPLowerInvisc, "--", color=cpLowColor, linewidth=1.0)
            axs[iAxsCP].set_ylim(self.CPlim)
            axs[iAxsCP].set_ylabel("$c_p$", rotation="horizontal", ha="right", va="center")

            # Make legend for viscous vs. inviscid
            customLines = [
                mpl_lines.Line2D([0], [0], linestyle="-", color="k"),
                mpl_lines.Line2D([0], [0], linestyle="--", color="k", linewidth=1.0),
            ]
            axs[iAxsCP].legend(customLines, ["Viscous", "Inviscid"])

            # Plot the skin friction coefficient
            axs[iAxsCF].plot([min(x) - 1, max(x) + 1], [0.0, 0.0], zorder=-2, color="k", linewidth=0.7, alpha=0.3)
            axs[iAxsCF].plot(xCFUpper, CFUpper, color="k", zorder=-1, alpha=0.15)
            axs[iAxsCF].plot(xCFLower, CFLower, color="k", zorder=-1, alpha=0.15)
            axs[iAxsCF].plot(xCFUpper, CFUpper, color=cpUpColor)
            axs[iAxsCF].plot(xCFLower, CFLower, color=cpLowColor)
            axs[iAxsCF].set_ylim(self.CFlim)
            axs[iAxsCF].set_ylabel("$c_f$", rotation="horizontal", ha="right", va="center")

            # Plot the airfoil on the lower axis
            axs[iAxsFoil].plot(x, y, color="k", zorder=-1, alpha=0.15)
            axs[iAxsFoil].plot(x, y, color=color)
            axs[iAxsFoil].set_xlim(self.xlimFoil)
            axs[iAxsFoil].set_ylim(self.ylimFoil)
            axs[iAxsFoil].set_xlabel("x/c")
            axs[iAxsFoil].set_ylabel("y/c", rotation="horizontal", ha="right", va="center")
            axs[iAxsFoil].set_aspect("equal")
            axs[iAxsFoil].spines["right"].set_visible(False)
            axs[iAxsFoil].spines["top"].set_visible(False)
            axs[iAxsFoil].yaxis.set_ticks_position("left")
            axs[iAxsFoil].xaxis.set_ticks_position("bottom")

            if fileName is None and showPlot:
                plt.pause(0.5)
            self.airfoilFig = fig
            self.airfoilAxs = axs
        else:
            self._updateAirfoilPlot(pause=showPlot)

        if fileName is not None:
            self.airfoilFig.savefig(fileName)

        return self.airfoilFig, self.airfoilAxs

    def _updateAirfoilPlot(self, pause=True):
        """
        Updates the airfoil plot with current airfoil shape.
        Assumes that the current figure is the one with the
        airfoil on it and it was the most recently plotted line.
        This function should not be called directly by the user
        unless you know what you're doing.

        Parameters
        ----------
        pause : bool
            If true, pauses after updating plot for 0.5 sec. This is necessary
            when updating a live plot, but breaks the animation in postprocess.
        """
        if plt is None:
            raise ImportError("matplotlib is required to use the plotting tools")

        # Get data to plot
        x = self.coords[:, 0]
        y = self.coords[:, 1]
        CPUpper = self.sliceData[self.curAP.name]["cp_visc_upper"]
        CPUpperInvisc = self.sliceData[self.curAP.name]["cp_invisc_upper"]
        xUpper = self.sliceData[self.curAP.name]["x_upper"]
        CPLower = self.sliceData[self.curAP.name]["cp_visc_lower"]
        CPLowerInvisc = self.sliceData[self.curAP.name]["cp_invisc_lower"]
        xLower = self.sliceData[self.curAP.name]["x_lower"]
        CFUpper = self.sliceData[self.curAP.name]["cf_upper"]
        xCFUpper = self.sliceData[self.curAP.name]["x_cf_upper"]
        CFLower = self.sliceData[self.curAP.name]["cf_lower"]
        xCFLower = self.sliceData[self.curAP.name]["x_cf_lower"]

        iAxsCP = 0
        iAxsCF = 1
        iAxsFoil = 2

        # Compute new plot limits
        stackedCP = np.hstack((CPUpper, CPUpperInvisc, CPLower, CPLowerInvisc))
        stackedCP = stackedCP[np.isfinite(stackedCP)]
        stackedCF = np.hstack((CFUpper, CFLower))
        stackedCF = stackedCF[np.isfinite(stackedCF)]

        CPlim = [max(stackedCP) + 0.2, min(stackedCP) - 0.2]
        CPlim = [max(CPlim[0], self.CPlim[0]), min(CPlim[1], self.CPlim[1])]  # inverted axis

        CFlim = [min(stackedCF) - 0.002, max(stackedCF) + 0.002]
        CFlim = [min(CFlim[0], self.CFlim[0]), max(CFlim[1], self.CFlim[1])]

        xlimFoil = [min(x) - 0.01, max(x) + 0.01]
        ylimFoil = [min(y) - 0.01, max(y) + 0.01]
        xlimFoil = [min(xlimFoil[0], self.xlimFoil[0]), max(xlimFoil[1], self.xlimFoil[1])]
        ylimFoil = [min(ylimFoil[0], self.ylimFoil[0]), max(ylimFoil[1], self.ylimFoil[1])]

        # CP plot
        self.airfoilAxs[iAxsCP].lines[-1].remove()
        self.airfoilAxs[iAxsCP].lines[-1].remove()
        self.airfoilAxs[iAxsCP].lines[-1].remove()
        self.airfoilAxs[iAxsCP].lines[-1].remove()
        self.airfoilAxs[iAxsCP].plot(xUpper, CPUpper, color=cpUpColor)
        self.airfoilAxs[iAxsCP].plot(xLower, CPLower, color=cpLowColor)
        self.airfoilAxs[iAxsCP].plot(xUpper, CPUpperInvisc, "--", color=cpUpColor, linewidth=1.0)
        self.airfoilAxs[iAxsCP].plot(xLower, CPLowerInvisc, "--", color=cpLowColor, linewidth=1.0)
        self.airfoilAxs[iAxsCP].set_ylim(CPlim)

        # CF plot
        self.airfoilAxs[iAxsCF].lines[-1].remove()
        self.airfoilAxs[iAxsCF].lines[-1].remove()
        self.airfoilAxs[iAxsCF].plot(xCFUpper, CFUpper, color=cpUpColor)
        self.airfoilAxs[iAxsCF].plot(xCFLower, CFLower, color=cpLowColor)
        self.airfoilAxs[iAxsCF].set_ylim(CFlim)

        self.airfoilAxs[iAxsFoil].lines[-1].remove()
        self.airfoilAxs[iAxsFoil].plot(x, y, color=color)
        self.airfoilAxs[iAxsFoil].set_xlim(xlimFoil)
        self.airfoilAxs[iAxsFoil].set_ylim(ylimFoil)

        if pause:
            plt.pause(1e-6)

    def _setCoordinates(self, x, y, setComplex=False):
        """
        Set the x and y coordinates in the compiled XFOIL layer.

        Parameters
        ----------
        x : ndarray (# pts,)
            x coordinates of airfoil
        y : ndarray (# pts,)
            y coordinates of airfoil
        setComplex : bool, optional
            Set coordinates of the complex version, otherwise will set
            the coordinates in the real version, by default False
        """
        N = 572  # This is VERY Important: The airfoil input must be a
        # FIXED length of 572. Simply set the coordinates up
        # to NB and leave the remainder as zeros
        NB = len(x)

        if setComplex:
            dtype = complex
            xfoil = self.xfoil_cs
        else:
            dtype = float
            xfoil = self.xfoil

        x_input = np.zeros(N, dtype=dtype)
        y_input = np.zeros(N, dtype=dtype)
        x_input[:NB] = np.array(x).copy()
        y_input[:NB] = np.array(y).copy()

        xfoil.ci04.nb = NB
        xfoil.cr14.xb = x_input
        xfoil.cr14.yb = y_input
        xfoil.xfoil()

    @staticmethod
    def _getDefaultOptions():
        return {
            "maxIters": [int, 100],  # maximum iterations for XFOIL solver
            "printRealConvergence": [bool, True],
            "printComplexConvergence": [bool, False],
            "writeCoordinates": [bool, True],  # whether to write coordinates to .dat file when `writeSolution` called
            "writeSliceFile": [bool, True],  # whether or not to save chordwise data
            "writeSolution": [bool, False],  # whether or not to call writeSolution after each call to XFOIL
            "plotAirfoil": [bool, False],  # whether to plot airfoil while running
            "outputDirectory": [str, "."],  # where to put output files
            "numberSolutions": [bool, True],  # whether to add call counter to output file names
            "xTrip": [
                np.ndarray,
                np.full(2, np.nan),
            ],  # boundary layer trip x coordinate of upper and lower surface, respectively (two-element array)
            "nCrit": [float, 9.0],
            "rhoKS": [float, 500.0],  # aggregation parameter
        }

    @staticmethod
    def _readDat(filename, headerlines=0):
        """
        Reads a '.dat' style airfoil coordinate file,
        with each coordinate on a new line and each
        line containing an xy pair separate by whitespace.

        Parameters
        ----------
        filename : str
            dat file from which to read data
        headerlines : int
            the number of lines to skip at the beginning of the file to reach the coordinates

        Returns
        -------
        X : Ndarray [N,2]
            The coordinates read from the file
        """
        with open(filename, "r") as f:
            for _i in range(headerlines):
                f.readline()
            r = []
            while True:
                line = f.readline()
                if not line:
                    break  # end of file
                if line.isspace():
                    break  # blank line
                r.append([float(s) for s in line.split()])

                X = np.array(r)

        return X

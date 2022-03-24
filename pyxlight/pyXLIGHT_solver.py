"""
pyXLIGHT

pyXLIGHT is a wrapper for Mark Drela's XFOIL code. The purpose of this
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

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
from baseclasses import BaseSolver
from prefoil.preFoil import readCoordFile, _writeDat
from pygeo.pyGeo import pyGeo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    import niceplots as nice

    nice.setRCParams()
    plt.rcParams["font.size"] = 18
    colors = nice.get_niceColors()
    color = colors["Blue"]
    cp_up_color = colors["Blue"]
    cp_low_color = colors["Red"]
except ImportError:
    print("Install niceplots for nice looking airfoil figures")
    color = "b"
    cp_up_color = "b"
    cp_low_color = "r"

# =============================================================================
# Extension modules
# =============================================================================
from . import MExt
from .pyXLIGHT import xfoilAnalysis


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
        self.DATFileName = fileName
        self.coords = readCoordFile(fileName)
        self.coords = np.hstack((self.coords, np.zeros((self.coords.shape[0], 1))))
        self.coords0 = self.coords.copy()  # initial coordinates (never changes)
        self.setCoordinates(self.coords0)  # set the initial coordinates
        self.pointSetKwargs = {}

        self.curAP = None
        self.DVGeo = None

        # Dictionary with dictionary of functions for each aero problem
        self.funcs = {}
        self.funcsComplex = {}
        self.functionList = ["cl", "cd", "cm"]  # available functions

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
        self.sliceData = {}
        self.sliceDataComplex = {}

        # Possible AeroProblem design variables (only alpha for pyXLIGHT)
        self.possibleAeroDVs = ["alpha"]

        # Figure and axes used by self.plotAirfoil and self.updateAirfoilPlot
        self.airfoilFig = None
        self.airfoilAxs = None

    def setDVGeo(self, DVGeo, pointSetKwargs=None):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that PYXLIGHT does not **strictly** need a
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
        ptSetName = "pyxlight_%s_coords" % aeroProblem.name

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
            aeroProblem.callCounter
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

    def setCoordinatesComplex(self, coords):
        """
        Update the complex airfoil coordinates and associated point sets.

        Parameters
        ----------
        coords : ndarray
            New airfoil coordinates (either 2 or 3 columns)
        """
        super().setCoordinatesComplex(coords[:, 0].astype(complex), coords[:, 1].astype(complex))

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
        xfoil = self.xfoil
        funcs = self.funcs
        sliceData = self.sliceData
        var_type = float
        if useComplex:
            var_type = complex
            xfoil = self.xfoil_cs
            funcs = self.funcsComplex
            sliceData = self.sliceDataComplex

        self.setAeroProblem(aeroProblem)

        # Set flight condition and options
        xfoil.cr15.reinf1 = aeroProblem.re  # Reynolds number
        xfoil.cr09.minf1 = aeroProblem.mach  # Mach Number set
        xfoil.cr09.adeg = aeroProblem.alpha
        xfoil.ci04.itmax = self.getOption("maxIters")  # Iterations Limit Set

        # Call XFOIL
        xfoil.oper()

        # Store results in dictionary for current aero problem
        funcs[aeroProblem.name] = {
            "cl": var_type(xfoil.cr09.cl),
            "cd": var_type(xfoil.cr09.cd),
            "cm": var_type(xfoil.cr09.cm),
        }

        # Pull out and process the pressure coefficient data
        cpv = xfoil.cr04.cpv  # viscous
        cpi = xfoil.cr04.cpi  # inviscid
        x = xfoil.cr05.x
        y = xfoil.cr05.y
        end_foil_idx = np.argmax(x > self.coords[0, 0]) + 1  # XFOIL includes wake panels, which we don't want
        idx_lower_start = end_foil_idx // 2  # first half of data is upper surface
        sliceData[aeroProblem.name] = {
            "cp_visc_upper": cpv[:idx_lower_start].copy().astype(var_type),
            "cp_invisc_upper": cpi[:idx_lower_start].copy().astype(var_type),
            "x_upper": x[:idx_lower_start].copy().astype(var_type),
            "y_upper": y[:idx_lower_start].copy().astype(var_type),
            "cp_visc_lower": cpv[idx_lower_start:end_foil_idx].copy().astype(var_type),
            "cp_invisc_lower": cpi[idx_lower_start:end_foil_idx].copy().astype(var_type),
            "x_lower": x[idx_lower_start:end_foil_idx].copy().astype(var_type),
            "y_lower": y[idx_lower_start:end_foil_idx].copy().astype(var_type),
        }

        # Check for failure
        self.curAP.solveFailed = self.curAP.fatalFail = xfoil.cl01.lexitflag != 0

        # If not a derivative call, increment callCounter
        if not deriv:
            self.curAP.callCounter += 1

        # Write solution files if desired
        if not deriv and self.getOption("writeSolution"):
            self.writeSolution()

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """This is a generic shell function that potentially writes the various output files. The intent is that the
        user or calling program can call this file and ADflow write all the files that the user has defined. It is
        recommended that this function is used along with the associated logical flags in  the options to determine the
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

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + "_%3.3d" % number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
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
        _writeDat(fileName, self.coords[:, 0], self.coords[:, 1])

    def writeSlice(self, fileName):
        """Write pickle file containing the sliceData dictionary. The data can be
        accessed using the AeroProblem name as the key. Within that is a dictionary containing:
            "cp_visc_upper": viscous CP on the airfoil's upper surface
            "cp_invisc_upper": inviscid CP on the airfoil's upper surface
            "x_upper": x coordinates of the upper surface CP data
            "y_upper": y coordinates of the upper surface CP data
            "cp_visc_lower": viscous CP on the airfoil's lower surface
            "cp_invisc_lower": inviscid CP on the airfoil's lower surface
            "x_lower": x coordinates of the lower surface CP data
            "y_lower": y coordinates of the lower surface CP data

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

    def checkAdjointFailure(self, aeroProblem, funcsSens):
        """
        Pass through to checkSolutionFailure to maintain the same interface as ADFLOW.

        Take in an aeroProblem and check for adjoint failure, Then append the
        fail flag in funcsSens. Information regarding whether or not the
        last analysis with the aeroProblem was sucessful is
        included. This information is included as "funcsSens['fail']". If
        the 'fail' entry already exits in the dictionary the following
        operation is performed:

        funcsSens['fail'] = funcsSens['fail'] or <did this problem fail>

        In other words, if any one problem fails, the funcsSens['fail']
        entry will be True. This information can then be used
        directly in multiPointSparse. For direct interface with pyOptSparse
        the fail flag needs to be returned separately from the funcs.

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
            if mode="CS")
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
        """This the main python gateway for producing forward mode jacobian
        vector products. It is not generally called by the user by
        rather internally or from another solver. A DVGeo object and a
        mesh object must both be set for this routine.
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
            if mode="CS")
        Returns
        -------
        dict
            Jacobian vector product of evalFuncs given perturbation
        """
        if xDvDot is None and xSDot is None:
            raise ValueError("xDvDot and xSDot cannot both be None")

        if mode not in ["FD", "CS"]:
            raise ValueError(f'Jacobian vector product mode "{mode}" invalid. Must be either "FD" or "CS"')

        geoDVs = list(self.DVGeo.getValues().keys())
        possibleDVs = self.possibleAeroDVs + geoDVs
        for DV in xDvDot.keys():
            if DV not in possibleDVs:
                raise ValueError(f'Perturbed design variable "{DV}" is not valid')

        if h is None:
            if mode == "FD":
                h = 1e-6
                orig_funcs = deepcopy(self.funcs[self.curAP.name])
            elif mode == "CS":
                h = 1e-100j

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

        Parameters
        ----------
        offsetDist : float
            Distance to extrude airfoil (same units as airfoil coordinates)

        Returns
        -------
        pyGeo surface
            Extruded airfoil surface
        """
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
            teHeight=0.25 * 0.0254,
        )

    def plotAirfoil(self, fileName=None):
        """
        Plots the current airfoil and returns the figure.

        Parameters
        ----------
        fileName : str, optional
            FileName to save to, if none specified it will
            show the plot with plt.show()

        Returns
        -------
        matplotlib figure
            Figure with airfoil plotted to it
        """

        if self.airfoilFig is None:
            # Get data to plot
            x = self.coords[:, 0]
            y = self.coords[:, 1]
            CPUpper = self.sliceData[self.curAP.name]["cp_visc_upper"]
            CPUpper_invisc = self.sliceData[self.curAP.name]["cp_invisc_upper"]
            xUpper = self.sliceData[self.curAP.name]["x_upper"]
            CPLower = self.sliceData[self.curAP.name]["cp_visc_lower"]
            CPLower_invisc = self.sliceData[self.curAP.name]["cp_invisc_lower"]
            xLower = self.sliceData[self.curAP.name]["x_lower"]

            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[10, 10])
            plt.ion()
            plt.show()

            # Plot the CP on the upper axis
            axs[0].plot(xUpper, CPUpper, color="k", zorder=-1, alpha=0.15)
            axs[0].plot(xLower, CPLower, color="k", zorder=-1, alpha=0.15)
            axs[0].plot(xUpper, CPUpper, color=cp_up_color)
            axs[0].plot(xLower, CPLower, color=cp_low_color)
            axs[0].plot(xUpper, CPUpper_invisc, "--", color=cp_up_color, linewidth=1.)
            axs[0].plot(xLower, CPLower_invisc, "--", color=cp_low_color, linewidth=1.)
            axs[0].invert_yaxis()
            axs[0].set_ylabel("$c_p$", rotation="horizontal", ha="right", va="center")

            # Make legend for viscous vs. inviscid
            customLines = [Line2D([0], [0], linestyle="-", color="k"),
                            Line2D([0], [0], linestyle="--", color="k", linewidth=1.)]
            axs[0].legend(customLines, ["Viscous", "Inviscid"])

            # Plot the airfoil on the lower axis
            axs[1].plot(x, y, color="k", zorder=-1, alpha=0.15)
            axs[1].plot(x, y, color=color)
            axs[1].set_xlim([min(x) - 0.01, max(x) + 0.01])
            axs[1].set_ylim([min(y) - 0.01, max(y) + 0.01])
            axs[1].set_xlabel("x/c")
            axs[1].set_ylabel("y/c", rotation="horizontal", ha="right", va="center")
            axs[1].set_aspect("equal")
            axs[1].spines["right"].set_visible(False)
            axs[1].spines["top"].set_visible(False)
            axs[1].yaxis.set_ticks_position("left")
            axs[1].xaxis.set_ticks_position("bottom")

            if fileName is None:
                plt.pause(0.5)
            self.airfoilFig = fig
            self.airfoilAxs = axs
        else:
            self.updateAirfoilPlot()

        if fileName is not None:
            self.airfoilFig.savefig(fileName)

        return self.airfoilFig

    def updateAirfoilPlot(self):
        """
        Updates the airfoil plot with current airfoil shape.
        Assumes that the current figure is the one with the
        airfoil on it and it was the most recently plotted line.
        """
        # Get data to plot
        y0 = self.coords0[:, 1]
        x = self.coords[:, 0]
        y = self.coords[:, 1]
        CPUpper = self.sliceData[self.curAP.name]["cp_visc_upper"]
        CPUpper_invisc = self.sliceData[self.curAP.name]["cp_invisc_upper"]
        xUpper = self.sliceData[self.curAP.name]["x_upper"]
        CPLower = self.sliceData[self.curAP.name]["cp_visc_lower"]
        CPLower_invisc = self.sliceData[self.curAP.name]["cp_invisc_lower"]
        xLower = self.sliceData[self.curAP.name]["x_lower"]

        # CP plot
        self.airfoilAxs[0].lines.pop(-1)
        self.airfoilAxs[0].lines.pop(-1)
        self.airfoilAxs[0].lines.pop(-1)
        self.airfoilAxs[0].lines.pop(-1)
        self.airfoilAxs[0].plot(xUpper, CPUpper, color=cp_up_color)
        self.airfoilAxs[0].plot(xLower, CPLower, color=cp_low_color)
        self.airfoilAxs[0].plot(xUpper, CPUpper_invisc, "--", color=cp_up_color, linewidth=1.)
        self.airfoilAxs[0].plot(xLower, CPLower_invisc, "--", color=cp_low_color, linewidth=1.)

        self.airfoilAxs[1].lines.pop(-1)
        self.airfoilAxs[1].plot(x, y, color=color)
        self.airfoilAxs[1].set_ylim([min(min(y), min(y0)) - 0.01, max(max(y), max(y0)) + 0.01])
        plt.pause(0.5)

    @staticmethod
    def _getDefaultOptions():
        return {
            "maxIters": [int, 100],  # maximum iterations for XFOIL solver
            "writeCoordinates": [bool, True],  # Whether to write coordinates to .dat file when `writeSolution` called
            "writeSliceFile": [bool, True],  # Whether or not to save chordwise data
            "writeSolution": [bool, False],  # Whether or not to call writeSolution after each call to XFOIL
            "plotAirfoil": [bool, False],  # Whether to plot airfoil while running
            "outputDirectory": [str, "."],  # Where to put output files
            "numberSolutions": [bool, True], # Whether to add call counter to output file names
        }

# ======================================================================
#         Import modules
# ======================================================================
# rst imports (beg)
import os
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVConstraints, DVGeometryCST
from pyoptsparse import Optimization, OPT
from multipoint import multiPointSparse
from cmplxfoil import CMPLXFOIL, AnimateAirfoilOpt

# rst imports (end)

# ======================================================================
#         Specify parameters for optimization
# ======================================================================
# rst params (beg)
mycl = 0.5  # lift coefficient constraint
alpha = 0.0 if mycl == 0.0 else 1.0  # initial angle of attack (zero if the target cl is zero)
mach = 0.1  # Mach number
Re = 1e6  # Reynolds number
T = 288.15  # 1976 US Standard Atmosphere temperature @ sea level (K)
# rst params (end)

# ======================================================================
#         Create multipoint communication object
# ======================================================================
# rst procs (beg)
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet("cruise", nMembers=1, memberSizes=MPI.COMM_WORLD.size)
MP.createCommunicators()
# rst procs (end)

# ======================================================================
#         Create output directory
# ======================================================================
# rst dir (beg)
curDir = os.path.abspath(os.path.dirname(__file__))
outputDir = os.path.join(curDir, "output")

if not os.path.exists(outputDir):
    os.mkdir(outputDir)
# rst dir (end)

# ======================================================================
#         CFD solver set-up
# ======================================================================
# rst solver (beg)
aeroOptions = {
    "writeSolution": True,
    "writeSliceFile": True,
    "writeCoordinates": True,
    "plotAirfoil": True,
    "outputDirectory": outputDir,
}

# Create solver
CFDSolver = CMPLXFOIL(os.path.join(curDir, "naca0012.dat"), options=aeroOptions)
# rst solver (end)

# ======================================================================
#         Set up flow conditions with AeroProblem
# ======================================================================
# rst ap (beg)
ap = AeroProblem(
    name="fc",
    alpha=alpha if mycl != 0.0 else 0.0,
    mach=mach,
    reynolds=Re,
    reynoldsLength=1.0,
    T=T,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=["cl", "cd"],
)
# Add angle of attack variable
if mycl != 0.0:
    ap.addDV("alpha", value=ap.alpha, lower=-10.0, upper=10.0, scale=1.0)
# rst ap (end)

# ======================================================================
#         Geometric Design Variable Set-up
# ======================================================================
# rst geom (beg)
DVGeo = DVGeometryCST(os.path.join(curDir, "naca0012.dat"))

nCoeff = 4  # number of CST coefficients on each surface
DVGeo.addDV("upper_shape", dvType="upper", dvNum=nCoeff, lower=-.1, upper=.5)
DVGeo.addDV("lower_shape", dvType="lower", dvNum=nCoeff, lower=-.5, upper=.1)

# Add DVGeo object to CFD solver
CFDSolver.setDVGeo(DVGeo)
# rst geom (end)

# ======================================================================
#         DVConstraint Setup
# ======================================================================
# rst cons (beg)
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Thickness, volume, and leading edge radius constraints
le = 0.0001
wingtipSpacing = 0.1
leList = [[le, 0, wingtipSpacing], [le, 0, 1.0 - wingtipSpacing]]
teList = [[1.0 - le, 0, wingtipSpacing], [1.0 - le, 0, 1.0 - wingtipSpacing]]
DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=0.85, scaled=True)
DVCon.addThicknessConstraints2D(leList, teList, 2, 100, lower=0.25, scaled=True)
le = 0.01
leList = [[le, 0, wingtipSpacing], [le, 0, 1.0 - wingtipSpacing]]
DVCon.addLERadiusConstraints(leList, 2, axis=[0, 1, 0], chordDir=[-1, 0, 0], lower=0.85, scaled=True)

fileName = os.path.join(outputDir, "constraints.dat")
DVCon.writeTecplot(fileName)
# rst cons (end)

# ======================================================================
#         Functions:
# ======================================================================
# rst funcs (beg)
def cruiseFuncs(x):
    print(x)
    # Set design vars
    DVGeo.setDesignVars(x)
    ap.setDesignVars(x)
    # Run CFD
    CFDSolver(ap)
    # Evaluate functions
    funcs = {}
    DVCon.evalFunctions(funcs)
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print("functions:")
        for key, val in funcs.items():
            if key == "DVCon1_thickness_constraints_0":
                continue
            print(f"    {key}: {val}")
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    CFDSolver.checkAdjointFailure(ap, funcsSens)
    print("function sensitivities:")
    evalFunc = ["fc_cd", "fc_cl", "fail"]
    for var in evalFunc:
        print(f"    {var}: {funcsSens[var]}")
    return funcsSens


def objCon(funcs, printOK):
    # Assemble the objective and any additional constraints:
    funcs["obj"] = funcs[ap["cd"]]
    funcs["cl_con_" + ap.name] = funcs[ap["cl"]] - mycl
    if printOK:
        print("funcs in obj:", funcs)
    return funcs


# rst funcs (end)

# ======================================================================
#         Optimization Problem Set-up
# ======================================================================
# rst optprob (beg)
# Create optimization problem
optProb = Optimization("opt", MP.obj)

# Add objective
optProb.addObj("obj", scale=1e4)

# Add variables from the AeroProblem
ap.addVariablesPyOpt(optProb)

# Add DVGeo variables
DVGeo.addVariablesPyOpt(optProb)

# Add constraints
DVCon.addConstraintsPyOpt(optProb)

# Add cl constraint
optProb.addCon("cl_con_" + ap.name, lower=0.0, upper=0.0, scale=1.0)

# Enforce first upper and lower CST coefficients to add to zero
# to maintain continuity at the leading edge
jac = np.zeros((1, nCoeff), dtype=float)
jac[0, 0] = 1.0
optProb.addCon(
    "first_cst_coeff_match",
    lower=0.0,
    upper=0.0,
    linear=True,
    wrt=["upper_shape", "lower_shape"],
    jac={"upper_shape": jac, "lower_shape": jac},
)

# The MP object needs the 'obj' and 'sens' function for each proc set,
# the optimization problem and what the objcon function is:
MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)
optProb.printSparsity()
optProb.getDVConIndex()
# rst optprob (end)

# rst opt (beg)
# Run optimization
optOptions = {"IFILE": os.path.join(outputDir, "SLSQP.out")}
opt = OPT("SLSQP", options=optOptions)
sol = opt(optProb, MP.sens, storeHistory=os.path.join(outputDir, "opt.hst"))
if MPI.COMM_WORLD.rank == 0:
    print(sol)
# rst opt (end)

# ======================================================================
#         Postprocessing
# ======================================================================
# rst postprocessing (beg)
# Save the final figure
CFDSolver.airfoilAxs[1].legend(["Original", "Optimized"], labelcolor="linecolor")
CFDSolver.airfoilFig.savefig(os.path.join(outputDir, "OptFoil.pdf"))

# Animate the optimization
AnimateAirfoilOpt(outputDir, "fc").animate(
    outputFileName=os.path.join(outputDir, "OptFoil"), fps=10, dpi=300, extra_args=["-vcodec", "libx264"]
)
# rst postprocessing (end)

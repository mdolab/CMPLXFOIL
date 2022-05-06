"""
This file contains tools for postprocessing airfoil optimizations
from the PYXLIGHT solver class in pyXLIGHT_solver.py.

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
import pickle as pkl
from glob import glob

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
from baseclasses import AeroProblem
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# =============================================================================
# Extension modules
# =============================================================================
from .pyXLIGHT_solver import PYXLIGHT


class AnimateAirfoilOpt:
    """
    Class for generating animations of airfoil optimization
    """

    def __init__(self, dirName, APName):
        """
        Initialize the object with the directory and AeroProblem name.
        This object assumes that the files are accessible under the name
        <dirName>/<APName>_<iteration number>.<dat or pkl>
        and that BOTH dat (airfoil shape) and pkl (chordwise data) files
        are available. It also assumes that the iteration number starts at 1.

        Parameters
        ----------
        dirName : str
            Name of directory that contains the airfoil optimization data files
        APName : str
            Name of the AeroProblem to be animated.
        """
        self.dirName = dirName
        self.APName = APName

        # Get the list of file names and figure out how many iterations there are to plot
        self.fileList = glob(os.path.join(self.dirName, f"{self.APName}_*.dat"))
        for i in range(len(self.fileList)):
            self.fileList[i] = self.fileList[i][:-4]
        self.fileList.sort()
        self.iters = len(self.fileList)

        # Make sure pickle files are also available
        if len(glob(os.path.join(self.dirName, f"{self.APName}_*.pkl"))) != self.iters:
            raise FileNotFoundError("There must be a pickle file and dat file for each iteration")

        print(f"Found dat and pkl files for {self.iters} iterations")

    def animate(self, outputFileName="airfoil_opt", ext="mp4", **animKwargs):
        """
        Generate an animation of an optimization.

        Parameters
        ----------
        outputFileName : str, optional
            Movie filename to save to with no extension (default "airfoil_opt")
        ext : str, optional
            Extension for animation ("mp4" and "gif" are useful ones)
        animKwargs : optional
            Additional keyword arguments to be passed to matplotlib's
            FuncAnimation save method
        """
        # If no animation keyword args specified and file type is mp4, set codec
        if not animKwargs:
            animKwargs = {"fps": 15, "dpi": 200}
            if ext.lower() == "mp4":
                animKwargs["extra_args"] = ["-vcodec", "libx264"]

        # Create initial plot
        foil = PYXLIGHT(self.fileList[0] + ".dat")
        foil.curAP = AeroProblem(self.APName, mach=0.5, altitude=0.0)
        with open(self.fileList[0] + ".pkl", "rb") as f:
            foil.sliceData = pkl.load(f)
        fig, axs = foil.plotAirfoil()
        CPlim = foil.CPlim
        CFlim = foil.CFlim
        xlimFoil = foil.xlimFoil
        ylimFoil = foil.ylimFoil
        coords0 = foil.coords0

        def animateFrame(i):
            """
            Function to be called by FuncAnimation
            """
            print(f"Rendering frame {i} of {self.iters}.....{(i + 1)/self.iters * 100:0.2f}% done", end="\r")
            foil = PYXLIGHT(self.fileList[i] + ".dat")
            foil.curAP = AeroProblem(self.APName, mach=0.5, altitude=0.0)
            foil.CPlim = CPlim
            foil.CFlim = CFlim
            foil.xlimFoil = xlimFoil
            foil.ylimFoil = ylimFoil
            foil.coords0 = coords0
            with open(self.fileList[i] + ".pkl", "rb") as f:
                foil.sliceData = pkl.load(f)
            foil.airfoilFig = fig
            foil.airfoilAxs = axs
            foil.updateAirfoilPlot(pause=False)

        # Call the animator and save the result as a movie file
        anim = FuncAnimation(fig, animateFrame, frames=self.iters, interval=66, blit=False)
        anim.save(outputFileName + "." + ext, **animKwargs)

        # Save the last frame
        animateFrame(self.iters - 1)
        plt.savefig(outputFileName + ".pdf")

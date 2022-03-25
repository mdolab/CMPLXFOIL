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
import time

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


class AnimateAirfoilOpt():
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

        # Figure out how many iterations there are to plot
        # This is a pretty bad algorithm for checking (N^2), but it is robust
        # to the way the iteration number is printed in the file and fast
        # compared to the animation.
        i = 0
        datExists = True
        pklExists = True
        while datExists and pklExists:
            i += 1
            datExists = False
            pklExists = False

            if self.findFile(i, 'dat') is not None:
                datExists = True
            
            if self.findFile(i, 'pkl') is not None:
                pklExists = True

        self.iters = i - 1
        print(f"Found dat and pkl files for {self.iters} iterations")
    
    def findFile(self, nIter, ext):
        """
        Find the file name associated with a given iteration.
        If the file is not found, returns None.

        Parameters
        ----------
        nIter : int
            Iteration number
        ext : str
            Extension (either "dat" or "pkl")
        """
        # See if the file is in the directory
        files = os.listdir(self.dirName)
        for f in files:
            # Get the iteration number and file type
            try:
                i = int(f.split('_')[-1].split('.')[0])
            except ValueError:
                continue
            fType = f.split('_')[-1].split('.')[1]
            AP = '_'.join(f.split('_')[:-1])

            if nIter == i and fType == ext and AP == self.APName:
                return os.path.join(self.dirName, f)

        return None

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
                animKwargs["extra_args"] = ['-vcodec', 'libx264']

        # Create initial plot
        foil = PYXLIGHT(self.findFile(1, 'dat'))
        foil.curAP = AeroProblem(self.APName, mach=0.5, altitude=0.)
        with open(self.findFile(1, 'pkl'), "rb") as f:
            foil.sliceData = pkl.load(f)
        fig, axs = foil.plotAirfoil()
        CPlim = foil.CPlim
        coords0 = foil.coords0

        def animateFrame(i):
            """
            Function to be called by FuncAnimation
            """
            print(f"Rendering frame {i} of {self.iters}.....{i/self.iters * 100:0.2f}% done", end="\r")
            i += 1  # the files are 1 indexed
            foil = PYXLIGHT(self.findFile(i, 'dat'))
            foil.curAP = AeroProblem(self.APName, mach=0.5, altitude=0.)
            foil.CPlim = CPlim
            foil.coords0 = coords0
            with open(self.findFile(i, 'pkl'), "rb") as f:
                foil.sliceData = pkl.load(f)
            foil.airfoilFig = fig
            foil.airfoilAxs = axs
            foil.updateAirfoilPlot()
            axs[0].invert_yaxis()
        
        # Call the animator and save the result as a movie file
        anim = FuncAnimation(fig, animateFrame,
                             frames=self.iters, interval=66, blit=False)
        anim.save(outputFileName + ".mp4", **animKwargs)

        # Save the last frame
        animateFrame(self.iters - 1)
        plt.savefig(outputFileName + ".pdf")

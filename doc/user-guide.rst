*********************
Optimization Tutorial
*********************

Introduction
============
This section describes a sample run script for airfoil optimization with CMPLXFOIL.
It is very similar to the `MACH-Aero single point airfoil tutorial <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machAeroTutorials/airfoilopt_singlepoint.html>`_.
This example uses pyOptSparse's SLSQP optimizer because it comes with pyOptSparse, but SNOPT is recommended for more robustness, speed, and tunability.

The optimization problem solved in this script is

| *minimize*
|    :math:`C_D`
| *with respect to*
|    4 upper surface shape variables (CST coefficients)
|    4 lower surface shape variables (CST coefficients)
| *subject to*
|    :math:`C_L = 0.5`
|    :math:`V \ge 0.85V_0`
|    :math:`R_{LE} \ge 0.85R_{LE, 0}`
|    :math:`t \ge 0.25t_0`
|    First upper surface CST coefficient = first lower surface CST coefficient

Getting set up
==============
In addition to the required CMPLXFOIL packages, this script requires `pygeo <https://github.com/mdolab/pygeo>`_, `multipoint <https://github.com/mdolab/multipoint>`_, `pyoptsparse <https://github.com/mdolab/pyoptsparse>`_, and mpi4py.
The script can be found in CMPLXFOIL/examples/single_point.py.

Dissecting the optimization script
==================================

Import libraries
----------------
.. literalinclude:: ../examples/single_point.py
    :start-after: # rst imports (beg)
    :end-before: # rst imports (end)

Specifying parameters for the optimization
------------------------------------------
These parameters define the flight condition and initial angle of attack for the optimization.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst params (beg)
    :end-before: # rst params (end)

Creating processor sets
-----------------------
Allocating sets of processors for different analyses can be helpful for multiple design points, but this is a single point optimization, so only one point is added.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst procs (beg)
    :end-before: # rst procs (end)

Creating output directory
-------------------------
This section creates a directory in the run script's directory in which to save files from the optimization.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst dir (beg)
    :end-before: # rst dir (end)

CMPLXFOIL solver setup
----------------------

The options tell the solver to write out chordwise aerodynamic data (``writeSliceFile``) and the airfoil coordinates (``writeCoordinates``) every time it is called.
It also enables live plotting during the optimization (``plotAirfoil``).
Finally, it specifies the output directory to save these files.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst solver (beg)
    :end-before: # rst solver (end)

Other options allow the user to adjust the maximum iterations allowed to the XFOIL solver and to specify a location at which to trip the boundary layer.

Set the AeroProblem
-------------------

We add angle of attack as a design variable (if the target lift coefficient is not zero) and set up the AeroProblem using given flow conditions.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst ap (beg)
    :end-before: # rst ap (end)

Geometric parametrization
-------------------------

This examples uses a class-shape transformation (CST) airfoil parameterization because it requires no additional files or other setup.
Four CST parameters are added to the upper and lower surface (the class shape and chord length are other possible design variables through DVGeometryCST).
The DVGeometryCST instance will set the initial design variables values by fitting them to the input dat file's geometry.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst geom (beg)
    :end-before: # rst geom (end)

Geometric constraints
---------------------
In this section, we add volume, thickness, and leading edge radius constraints.
They are chosen to achieve practical airfoils and to guide the optimizer away from infeasible design, such as the upper and lower surfaces crossing over each other.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst cons (beg)
    :end-before: # rst cons (end)

Optimization callback functions
-------------------------------
This section defines callback functions that are used by the optimizer to get objective, constraint, and derivative information.
See the `MACH-Aero <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machAeroTutorials/opt_aero.html>`_ aerodynamic optimization tutorial for more information.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst funcs (beg)
    :end-before: # rst funcs (end)

Optimization problem
--------------------
This section sets up the optimization problem by adding the necessary design variables and constraints.
An additional constraint for this problem enforces that the first upper and lower surface CST coefficients are equal.
This is to maintain continuity on the leading edge.
It also prints out some useful information about the optimization problem setup.
See the `MACH-Aero <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machAeroTutorials/opt_aero.html>`_ aerodynamic optimization tutorial for more information.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst optprob (beg)
    :end-before: # rst optprob (end)

Run optimization
----------------
Run the optimization using pyOptSparse's SLSQP optimizer and print the solution.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst opt (beg)
    :end-before: # rst opt (end)

Postprocessing
--------------
Finally, we save the final figure and use the built-in animation utility to create an optimization movie.

.. literalinclude:: ../examples/single_point.py
    :start-after: # rst postprocessing (beg)
    :end-before: # rst postprocessing (end)

Run it yourself!
================

To run the script, use the following command:

.. prompt:: bash

    python single_point.py

In the output directory, it should create the following animation after the optimization completes:

.. image:: assets/example_anim.gif
    :width: 600
    :align: center

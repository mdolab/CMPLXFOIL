pyXLIGHT
========
[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status/mdolab.pyXLIGHT?repoName=mdolab%2FpyXLIGHT&branchName=main)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=40&repoName=mdolab%2FpyXLIGHT&branchName=main)
[![Documentation Status](https://readthedocs.com/projects/mdolab-pyxlight/badge/?version=latest&token=7a9e7987d2288b741e09686619f4cd425b1a7348ebbcca59c0d20b2ad5a003f6)](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/?badge=latest)

pyXLIGHT is a version of Mark Drela's XFOIL code with the GUI features removed.
Gradient computation is implemented with the complex-step method.

Documentation
-------------
Please see the [documentation](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/) for installation and usage details.

To locally build the documentation, enter the ``doc`` folder and run ``make html`` in the command line prompt.
You can then view the built documentation in the ``doc/_build/html/`` by opening ``index.html``.

Citing pyXLIGHT
---------------
If you use pyXLIGHT, please see [this page](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/citation.html) for citation information.

License
-------
Copyright 2021 MDO Lab

Distributed using the GNU General Public License (GPL), version 2.0; see the LICENSE.md file for details.

## Methods to implement for pyXLIGHT `AeroSolver` class
- [ ] `addSlices`: Not exactly sure what this one is needed for
- [x] `setDVGeo`: Add `DVGeo` object to class
- [x] `setAeroProblem`: Do that
- [x] `setCoordinates`: yaaa
- [x] `getCoordinates`: also that
- [x] `getTriangulatedMeshSurface`: Return triangulated surface for `DVCon` use
- [x] `__call__`: Take in aero problem and run the solver
- [x] `evalFunctions`: Evaluate desired functions
- [x] `checkSolutionFailure`: Check for solution failure so optimizer knows about solver convergence
- [x] `evalFunctionSens`: Evaluate gradients of desired functions
- [x] `computeJacobianVectorProductFwd`
- [x] `checkAdjointFailure`: Should we return a failure if pyXLIGHT doesn't converge one of the complex step evaluations?

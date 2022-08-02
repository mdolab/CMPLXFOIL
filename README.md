<h2 align="center">
    <img src="/doc/assets/cmplxfoil_logo.svg" width="400" />
</h2>

[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status/mdolab.CMPLXFOIL?repoName=mdolab%2FCMPLXFOIL&branchName=main)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=45&repoName=mdolab%2FCMPLXFOIL&branchName=main)
[![Documentation Status](https://readthedocs.com/projects/mdolab-pyxlight/badge/?version=latest&token=7a9e7987d2288b741e09686619f4cd425b1a7348ebbcca59c0d20b2ad5a003f6)](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/?badge=latest)

CMPLXFOIL is a version of Mark Drela's XFOIL code with the GUI features removed.
Gradient computation is implemented with the complex-step method.

<p align="center">
  <img src="/doc/assets/airfoil_opt.gif" width="500">
</p>

Documentation
-------------
Please see the [documentation](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/) for installation and usage details.

To locally build the documentation, enter the ``doc`` folder and run ``make html`` in the command line prompt.
You can then view the built documentation in the ``doc/_build/html/`` by opening ``index.html``.

Citing CMPLXFOIL
---------------
If you use CMPLXFOIL, please see [this page](https://mdolab-pyxlight.readthedocs-hosted.com/en/latest/citation.html) for citation information.

License
-------
Copyright 2021 MDO Lab

Distributed using the GNU General Public License (GPL), version 2.0; see the LICENSE.md file for details.

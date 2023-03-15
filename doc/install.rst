************
Installation
************

Conda Installation
==================
If you use Anaconda/Conda, the simplest way to install CMPLXFOIL is through the `conda forge <https://anaconda.org/conda-forge/cmplxfoil>`_:

.. prompt:: bash

    conda install -c conda-forge cmplxfoil

Installation from Source
========================
Building and installing CMPLXFOIL from source is a multi-step process that is machine specific.
This guide goes through the required steps to compile and install CMPLXFOIL.


Requirements
------------
To compile the required XFOIL core components, a Fortran and a C compiler must be installed on your system.
This can be GNU / GFortran or Intel, default configuration files are packaged for both types of compilers.

In addition to standard compilers, CMPLXFOIL requires the following dependencies:

=================== ======= =======
Package             Version Notes
=================== ======= =======
Python              3.X.X
NumPy               ---     ``conda install numpy`` or ``pip install numpy``
`baseclasses`_      ---     ``pip install mdolab-baseclasses``
`pyGeo`_            ---     optional; required for ``getTriangulatedMeshSurface`` method
matplotlib          ---     optional; required for ``plotAirfoil`` method (``pip install matplotlib``)
niceplots           >=2.0.0 optional; recommended for ``plotAirfoil`` method (``pip install niceplots``)
=================== ======= =======

.. _baseclasses: https://github.com/mdolab/baseclasses
.. _pyGeo: https://github.com/mdolab/pygeo

Build and Installation
----------------------
Building CMPLXFOIL is handeled automatically by a set of Makefiles which are distributed with the code.
These Makefiles require configuration files which specify machine-specific parameters, such as compiler locations and flags.
Default configuration files for Linux GCC and Linux Intel are included in the ``config/defaults`` directory.
Copy a configuration file to the main ``config/`` folder using the command below and modify its contents for your system and installation.

.. prompt:: bash

    cp config/defaults/config.<version>.mk config/config.mk

Once the configuration file is adjusted as needed, CMPLXFOIL can be built by running ``make`` in the root directory:

.. prompt:: bash

    make

This will compile both the real and complex versions of CMPLXFOIL, generating Python libraries which reference the XFOIL Fortran modules.
These will be automatically copied to the ``cmplxfoil/`` directory.

Once the Python libraries are generated, install CMPLXFOIL by running pip install in the root directory:

.. prompt:: bash

    pip install .

Verification
------------
Tests are located in the ``tests/`` directory and can be run with the command:

.. prompt:: bash

    testflo -v .

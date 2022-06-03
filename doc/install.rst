Installation
============
CMPLXFOIL is a Python wrapper for the Fortran-based XFOIL code.
Building and installing the code is a multi-step process that is machine specific.
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
=================== ======= =======

To use the :ref:`CMPLXFOIL<CMPLXFOIL API>` interface for the `MACH-Aero <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/index.html>`_ framework, the following additional dependencies are required:

=================== ======= =======
Package             Version Notes
=================== ======= =======
`baseclasses`_      ---     ``pip install mdolab-baseclasses``
`preFoil`_          ---
`pyGeo`_            --- 
=================== ======= =======

.. _baseclasses: https://github.com/mdolab/baseclasses
.. _preFoil: https://github.com/mdolab/prefoil
.. _pyGeo: https://github.com/mdolab/pygeo

Finally, CMPLXFOIL's plotting features use

=================== ======= =======
Package             Version Notes
=================== ======= =======
matplotlib          ---     ``pip install matplotlib``
`niceplots`_        ---     optional
=================== ======= =======

.. _niceplots: https://github.com/mdolab/niceplots

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

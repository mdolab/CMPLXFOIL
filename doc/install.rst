Installation
============
pyXLIGHT is a Python wrapper for the Fortran-based XFOIL code.
Building and installing the code is a multi-step process that is machine specific.
This guide goes through the required steps to compile and install pyXLIGHT.


Requirements
------------
To compile the required XFOIL core components, a Fortran and a C compiler must be installed on your system.
This can be GNU / GFortran or Intel, default configuration files are packaged for both types of compilers.

In addition to standard compilers, pyXLIGHT requires the following dependencies:

=================== =======
Package             Version
=================== =======
Python              3.X.X
NumPy               ---
Sphinx              ---
sphinx-mdolab-theme ---
=================== =======

With the exception of Python, which must be installed on your system, these packages can be installed using:

.. code-block:: bash

    $ pip install numpy sphinx sphinx-mdolab-theme

Build and Installation
----------------------
Building pyXLIGHT is handeled automatically by a set of Makefiles which are distributed with the code.
These Makefiles require configuration files which specify machine-specific parameters, such as compiler locations and flags.
Default configuration files for Linux GCC and Linux Intel are included in the ``config/defaults`` directory.
Copy a configuration file to the main ``config/`` folder using the command below and modify its contents for your system and installation.

.. code-block:: bash

    $ cp config/defaults/config.<version>.mk config/config.mk

Once the configuration file is adjusted as needed, pyXLIGHT can be built by running ``make`` in the root directory:

.. code-block:: bash

    $ make

This will compile both the real and complex versions of pyXLIGHT, generating Python libraries which reference the XFOIL Fortran modules.
These will be automatically copied to the ``pyxlight/`` directory.

Once the Python libraries are generated, install pyXLIGHT by running pip install in the root directory:

.. code-block:: bash

    $ pip install .

Verification
------------
pyXLIGHT has a single functionality test, which can be run to ensure that the code can be run and executed.
This test is located in the ``tests/`` directory and can be run with the command:

.. code-block:: bash

    $ testflo -v .

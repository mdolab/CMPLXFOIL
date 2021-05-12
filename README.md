pyXLIGHT
========
[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status/mdolab.pyXLIGHT?repoName=mdolab%2FpyXLIGHT&branchName=master)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=40&repoName=mdolab%2FpyXLIGHT&branchName=master)

pyXLIGHT is a version of Mark Drela's XFOIL code with the GUI features removed.
Gradient computation is implemented with the complex-step method.

Installation
------------
To install pyXLIGHT, the XFOIL components must be compiled and wrapped with f2py.
The package has configuration files for Linux using GFortran and Intel compilers, located in the `config/defaults/` directory.
Copy one of the defaults to the base `config/` directory and adjust it as needed:

```
cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk
```

Then, build the Fortran code and f2py wrapper using the provided Makefile.
In the root directory, use the `make` command:
```
make
```

This will build both the real and complex version of the code, copying one Python library for each to the `pyxlight/` directory.
The Python package can then be installed as:

```
pip install .
```

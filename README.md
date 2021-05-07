pyXLIGHT
========

pyXLIGHT is a version of Mark Drela's XFOIL code with the GUI features removed.
Gradient computation is implemented with the complex-step method.

Installation
------------
To install pyXLIGHT, the XFOIL components must be compiled first.
The package has configuration files for GFortran and Intel compilers, which can be used to build the Fortran components:

Intel:
```
make intel
```

GFortran:
```
make gfortran
```

This will build both the real and complex version of the code, copying one Python library for each to the `pyxlight/` directory.
The Python package can then be installed as:

```
pip install .
```

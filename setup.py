from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("pyxlight/__init__.py").read(),
)[0]

setup(
    name="pyXLIGHT",
    version=__version__,
    description="A Python wrapped version of Mark Drela's XFOIL code with the GUI features removed.",
    long_description="""
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
      """,
    long_description_content_type="text/markdown",
    keywords="",
    author="",
    author_email="",
    url="https://github.com/mdolab/pyXLIGHT",
    license="",
    packages=[
        "pyxlight",
    ],
    package_data={"pyxlight": ["*.so"]},
    install_requires=[
        "numpy",
    ],
    extras_require={"testing": ["testflo"]},
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)

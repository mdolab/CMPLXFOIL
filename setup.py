from setuptools import setup
import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("pyxlight/__init__.py").read(),
)[0]

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pyXLIGHT",
    version=__version__,
    description="A Python wrapped version of Mark Drela's XFOIL code with the GUI features removed.",
    long_description=long_description,
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
    install_requires=["numpy", "mdolab-baseclasses", "prefoil"],
    extras_require={"testing": ["testflo"], "docs": ["sphinx-mdolab-theme"], "plotting": ["matplotlib"]},
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)

# Config file for Linux and GFortran

AR       = ar
AR_FLAGS = -rvs

FF90       = gfortran
FF90_FLAGS = -O2 -fPIC -fdefault-real-8 -fdefault-double-8

CC = gcc
CC_FLAGS = -O2 -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS = -fPIC -shared

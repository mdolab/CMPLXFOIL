# Config file for Linux and Intel compilers

AR       = ar
AR_FLAGS = -rvs

# Note that ";" is there to avoid make shell optimization, otherwise the shell command may fail
ICC_EXISTS := $(shell command -v icc;)

ifdef ICC_EXISTS
  # icc only exists on older Intel versions
  # Assume that we want to use the old compilers
  FF90       = ifort
  CC = icc
  OPT_FLAGS = -O2
else
  # Use the new compilers
  FF90       = ifx
  CC = icx
  OPT_FLAGS = -O0
endif

FF90_FLAGS = $(OPT_FLAGS) -r8 -fPIC
CC_FLAGS = $(OPT_FLAGS) -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS = -fPIC -shared

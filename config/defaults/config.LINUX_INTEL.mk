# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
# Note that ";" is there to avoid make shell optimization, otherwise the shell command may fail
ICC_EXISTS := $(shell command -v icc;)
ifdef ICC_EXISTS
  FF90       = ifort
else
  FF90       = ifx
endif

FF90_FLAGS = -O2 -r8 -fPIC

F2PY = f2py
F2PY_FF90 = intelem

# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
# Note that ";" is there to avoid make shell optimization, otherwise the shell command may fail
IFORT_EXISTS := $(shell command -v ifort;)
ifdef IFORT_EXISTS
  FF90       = ifort
  OPT_FLAGS  = -O2
else
  FF90       = ifx
  # We see a lot of solver failures with ifx when using anything more aggressive than -O0 optimization, increase the optimization level at your own risk
  OPT_FLAGS  = -O0
endif

FF90_FLAGS = $(OPT_FLAGS) -r8 -fPIC

F2PY = f2py
F2PY_FF90 = intelem

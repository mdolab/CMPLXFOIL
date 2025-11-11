# Config File for LINUX and GFortran Compiler
AR       = ar
AR_FLAGS = -rvs

FF90       = gfortran
FF90_FLAGS = -O2 -fdefault-real-8 -fPIC

F2PY = f2py
# Note: The F2PY_FF90 variable is ignored for numpy>=2.0 or python>=3.12. FF90 is used instead.
F2PY_FF90 = gnu95

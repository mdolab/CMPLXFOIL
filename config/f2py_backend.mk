# The following checks the numpy version and uses it to support both legacy distutils and meson backends
# Also python 3.12 and later removed distutils so we need to check for that as well
PYTHON_NEW_ENOUGH := $(shell python -c "import sys; print(int(sys.version_info[:2] >= (3, 12)))")
NUMPY_NEW_ENOUGH := $(shell python -c "import numpy; print(int(int(numpy.__version__.split('.')[0]) >= 2))")

# If either is true, use Meson
USE_MESON := 0
ifeq ($(PYTHON_NEW_ENOUGH),1)
    USE_MESON := 1
endif
ifeq ($(NUMPY_NEW_ENOUGH),1)
	USE_MESON := 1
endif

ifeq ($(USE_MESON),1)
	F2PY_CMD = export FC=$(FF90); $(F2PY)
	REAL_F2PY_LINK = -L$(shell pwd) -lxfoil
	COMPLEX_F2PY_LINK = -L$(shell pwd) -lxfoil_cs
else
	F2PY_CMD = $(F2PY) --fcompiler=$(F2PY_FF90)
	REAL_F2PY_LINK = libxfoil.a
	COMPLEX_F2PY_LINK = libxfoil_cs.a
endif

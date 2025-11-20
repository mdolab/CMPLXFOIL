# The following checks the numpy version and uses it to support both legacy distutils and meson backends
# Also python 3.12 and later removed distutils so we need to check for that as well
USE_MESON := $(shell python -c "import sys, numpy; print(int(sys.version_info[:2] >= (3,12) or int(numpy.__version__.split('.')[0]) >= 2))")

ifeq ($(USE_MESON),1)
	F2PY_CMD = export FC=$(FF90); $(F2PY)
	REAL_F2PY_LINK = -L$(shell pwd) -lxfoil
	COMPLEX_F2PY_LINK = -L$(shell pwd) -lxfoil_cs
else
	F2PY_CMD = $(F2PY) --fcompiler=$(F2PY_FF90)
	REAL_F2PY_LINK = libxfoil.a
	COMPLEX_F2PY_LINK = libxfoil_cs.a
endif

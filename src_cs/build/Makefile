# Include user supplied configuration file
include ../../config/config.mk

# Group Fortran compiler flags
FF90_ALL_FLAGS = $(FF90_FLAGS)

# Include full list of files
include fileList

# Include full list of directories
include directoryList

# Include the list of rules
include rules

# Set the make VPATH variable to the "dirs" variable from
# directorylist. We must first append the '../'
dirs:=$(addprefix ../,$(dirs))
VPATH:=$(dirs)

# Add ../ prefix to all source files, remove directories
fortranFiles:=$(addprefix ../,$(fortranFiles))
f77Files:=$(addprefix ../,$(f77Files))

fortranFilesNoDir=$(notdir $(fortranFiles))
f77FilesNoDir=$(notdir $(f77Files))

# Generate two separate list of .F90 and .f90 files using the filter command
f90Files=$(filter %.f90,$(fortranFilesNoDir))
F90Files=$(filter %.F90,$(fortranFilesNoDir))

# Finally convert all source files to .o
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(f77FilesNoDir:%.f=%.o)

.PHONY: sources lib

# The following checks the numpy version and uses it to support both legacy distutils and meson backends
# Also python 3.12 and later removed distutils so we need to check for that as well
USE_MESON := $(shell python -c "import sys, numpy; print(int(sys.version_info[:2] >= (3,12) or int(numpy.__version__.split('.')[0]) >= 2))")
F2PY_COMMON := --f90flags='$(FF90_ALL_FLAGS)' -c -m libcmplxfoil_cs ../f2py/libcmplxfoil_cs.pyf  -L$(shell pwd) -lxfoil_cs

default: lib ../f2py/libcmplxfoil_cs.pyf
ifeq ($(USE_MESON),1)
	export FC=$(FF90); \
	$(F2PY) $(F2PY_COMMON)
else
	$(F2PY) --fcompiler=$(F2PY_FF90) $(F2PY_COMMON)
endif
	mv libcmplxfoil_cs.*so ../../cmplxfoil/libcmplxfoil_cs.so

sources: $(OFILES)

lib: sources
	$(AR) $(AR_FLAGS) libxfoil_cs.a $(OFILES)

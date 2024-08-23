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

LIBDIR = .

LINKER_ALL_FLAGS = $(FF90_FLAGS) -L$(LIBDIR) -lcmplxfoil_cs $(LINKER_FLAGS)

PYTHON_OBJECTS = fortranobject.o \
		 libcmplxfoil_csmodule.o\
		 libcmplxfoil_cs-f2pywrappers.o\

default: lib ../f2py/libcmplxfoil_cs.pyf

# Generate Python inlude directory
	$(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())'))
	 @echo "#------------------------------------------------------#"
	 @echo Numpy Include Directory: $(NUMPY_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) ../../src/f2py/get_f2py.py))
	 @echo "#------------------------------------------------------#"
	 @echo f2py root directory: $(F2PY_ROOT)
	 @echo "#------------------------------------------------------#"

# Run f2py
	$(F2PY) ../f2py/libcmplxfoil_cs.pyf

# Compile f2py-generated C wrapper file
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libcmplxfoil_csmodule.c

# Compile the Fortran wrapper in the f2py source
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) -c \
	$(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compile f2py-generated Fortran wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c libcmplxfoil_cs-f2pywrappers.f

# Create shared library
	$(FF90) $(SO_LINKER_FLAGS) $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o libxfoil_cs.so
	mv libxfoil_cs.so libcmplxfoil_cs.so

# Test if the shared library can be imported in Python
	$(PYTHON) ../../src/f2py/importTest.py libcmplxfoil_cs.so

# Move the shared library to the main Python directory
	mv libcmplxfoil_cs.so ../../cmplxfoil

# Clean up temporary files
	-rm *.o
	-rm *.c
	-rm *.f

sources: $(OFILES)

lib: sources
	$(AR) $(AR_FLAGS) libcmplxfoil_cs.a $(OFILES)

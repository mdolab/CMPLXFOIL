#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************
# Include user supplied configuration file
include ../../config/config.mk

# Group Fortran compiler flags
FF_ALL_FLAGS = $(FFLAGS)

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
f77Files:=$(addprefix ../,$(f77Files))
f77FilesNoDir=$(notdir $(f77Files))

# Finally convert all source files to .o
OFILES=$(f77FilesNoDir:%.f=%.o)

.PHONY: sources lib

default: lib ../f2py/libpyxlight.pyf
	f2py  --fcompiler=$(FF90) --f90flags='$(FFLAGS)' -c -m libpyxlight ../f2py/libpyxlight.pyf libxfoil.a

sources: $(OFILES)

lib: sources
	$(AR) $(AR_FLAGS) libxfoil.a $(OFILES)

#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************
INCLUDE_COMMON = ../common.mk
include ${INCLUDE_COMMON}

XFOIL_LIB = libxfoil_cs.a
OBJ =	complexify.o\
	c_xfoil.o\
	c_xpanel.o\
	c_xoper.o\
	c_xsolve.o\
	c_xgdes.o\
	c_xbl.o\
	c_xblsys.o\
	c_xgeom.o\
	c_xutils.o\
	c_aread.o\
	c_userio.o\
	c_spline.o\
	c_naca.o

default: libxfoil_cs

clean:
	$(RM) *.o libxfoil_cs.a *.mod

libxfoil_cs: $(OBJ) 
	$(AR) $(AR_FLAGS) $(XFOIL_LIB) $(OBJ)


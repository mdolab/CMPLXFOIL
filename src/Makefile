#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************
INCLUDE_COMMON = ../common.mk
include ${INCLUDE_COMMON}

XFOIL_LIB = libxfoil.a
OBJ =	xfoil.o\
	 xpanel.o\
	 xoper.o\
	 xsolve.o\
	 xgdes.o\
	 xbl.o\
	 xblsys.o\
	 xgeom.o\
	 xutils.o\
	 aread.o\
	 userio.o\
	 spline.o\
	 naca.o\
	 xdriver.o\

default: libxfoil

clean:
	$(RM) *.o libxfoil.a

libxfoil: $(OBJ) 
	$(AR) $(AR_FLAGS) $(XFOIL_LIB) $(OBJ)


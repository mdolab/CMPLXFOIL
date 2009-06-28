#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************

#*********************************************************

##--------------------------
### Intel Fortran Compiler
#FC = ifort
#FFLAGS = -O2 -r8

##-------------------------
### Gfortran Compiler
FC = gfortran
FFLAGS = -O2 --default-real-8

SHELL = sh
BINDIR = .

AR       = ar
AR_FLAGS = -rvs
XFOIL_LIB = libxfoil_cs.a
PROGS = xfoil

SRC = .

XFOILOBJ =complexify.o c_xfoil.o  c_xpanel.o c_xoper.o c_xsolve.o c_xgdes.o c_xbl.o c_xblsys.o c_xgeom.o c_xutils.o c_aread.o c_userio.o c_spline.o c_naca.o
FTNLIB = 

default: libxfoil_cs

libxfoil_cs: ofiles
	 @echo " Adding to library ..."
	 @echo $(XFOIL_LIB)
	 @echo
	 $(AR) $(AR_FLAGS) $(XFOIL_LIB) \
	 $(XFOILOBJ) $(FNTLIB)
	 @echo 
	 @echo "library $(XFOIL_LIB) created."
	 @echo

clean:
	-/bin/rm *.o
	-/bin/rm libxfoil_cs.a


ofiles: $(XFOILOBJ)

complexify.o: $(SRC)/complexify.f90
	$(FC) -c $(FFLAGS)  $(SRC)/complexify.f90

c_xfoil.o: $(SRC)/c_xfoil.f $(SRC)/c_XFOIL.INC 
	$(FC) -c $(FFLAGS) $(SRC)/c_xfoil.f

c_xpanel.o: $(SRC)/c_xpanel.f  $(SRC)/c_XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xpanel.f

c_xoper.o: $(SRC)/c_xoper.f  $(SRC)/c_XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xoper.f

c_xsolve.o: $(SRC)/c_xsolve.f  $(SRC)/c_XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xsolve.f

c_xgdes.o: $(SRC)/c_xgdes.f  $(SRC)/c_XFOIL.INC 
	$(FC) -c $(FFLAGS)  $(SRC)/c_xgdes.f

c_xbl.o:   $(SRC)/c_xbl.f    $(SRC)/c_XFOIL.INC $(SRC)/c_XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xbl.f

c_xblsys.o: $(SRC)/c_xblsys.f                   $(SRC)/c_XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xblsys.f

c_xgeom.o: $(SRC)/c_xgeom.f $(SRC)/c_XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/c_xgeom.f

c_xutils.o: $(SRC)/c_xutils.f 
	$(FC) -c $(FFLAGS)  $(SRC)/c_xutils.f

c_aread.o: $(SRC)/c_aread.f 
	$(FC) -c $(FFLAGS)  $(SRC)/c_aread.f

c_userio.o: $(SRC)/c_userio.f 
	$(FC) -c $(FFLAGS)  $(SRC)/c_userio.f

c_spline.o: $(SRC)/c_spline.f 
	$(FC) -c $(FFLAGS)  $(SRC)/c_spline.f

c_naca.o: $(SRC)/c_naca.f 
	$(FC) -c $(FFLAGS)  $(SRC)/c_naca.f
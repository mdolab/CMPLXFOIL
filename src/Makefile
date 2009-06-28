#*********************************************************
# Makefile for XFOIL V6.93 programs
# H.Youngren 4/24/01
# M.Drela
#*********************************************************

##--------------------------
### Intel Fortran Compiler
#FC = ifort
#FFLAGS = -O3 -r8 -align all

##-------------------------
### Gfortran Compiler
FC = gfortran
FFLAGS = -O2 --default-real-8

SHELL = sh
BINDIR = .

AR       = ar
AR_FLAGS = -rvs
XFOIL_LIB = libxfoil.a
PROGS = xfoil

SRC = .

XFOILOBJ = xfoil.o xpanel.o xoper.o xsolve.o xgdes.o xbl.o xblsys.o xgeom.o xutils.o aread.o userio.o spline.o naca.o


default: libxfoil

libxfoil: ofiles
	 @echo " Adding to library ..."
	 @echo $(XFOIL_LIB)
	 @echo
	 $(AR) $(AR_FLAGS) $(XFOIL_LIB) \
	 $(XFOILOBJ) 
	 @echo 
	 @echo "library $(XFOIL_LIB) created."
	 @echo

clean:
	-/bin/rm *.o
	rm libxfoil.a

ofiles: $(XFOILOBJ)

xfoil.o: $(SRC)/xfoil.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xfoil.f

xpanel.o: $(SRC)/xpanel.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xpanel.f

xoper.o: $(SRC)/xoper.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xoper.f

xsolve.o: $(SRC)/xsolve.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xsolve.f

xgdes.o: $(SRC)/xgdes.f  $(SRC)/XFOIL.INC 
	$(FC) -c $(FFLAGS)  $(SRC)/xgdes.f

xbl.o:   $(SRC)/xbl.f    $(SRC)/XFOIL.INC $(SRC)/XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xbl.f

xblsys.o: $(SRC)/xblsys.f                  $(SRC)/XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xblsys.f

xgeom.o: $(SRC)/xgeom.f $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xgeom.f

xutils.o: $(SRC)/xutils.f 
	$(FC) -c $(FFLAGS)  $(SRC)/xutils.f

aread.o: $(SRC)/aread.f 
	$(FC) -c $(FFLAGS)  $(SRC)/aread.f

userio.o: $(SRC)/userio.f 
	$(FC) -c $(FFLAGS)  $(SRC)/userio.f

spline.o: $(SRC)/spline.f 
	$(FC) -c $(FFLAGS)  $(SRC)/spline.f

naca.o: $(SRC)/naca.f 
	$(FC) -c $(FFLAGS)  $(SRC)/naca.f

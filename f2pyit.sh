#/bin/bash
#Make script for pyXLIGHT

#/bin/bash
#Make script for pyXLIGHT

#First make the REAL source files
cd src

make
if [ ! $? -eq 0 ]; then
    exit
fi
cd ../

#Now f2py it!

# Intel Version
#f2py  --fcompiler=intel    -c -m pyxlight src/pyxlight.pyf src/libxfoil.a 

# gfortran Version
f2py  --fcompiler=gfortran -c -m pyxlight src/pyxlight.pyf src/libxfoil.a 


#Test the Module
python module_import_test.py real

sleep 2

#Next make the COMPLEX source files
cd src_cs

make
if [ ! $? -eq 0 ]; then
    exit
fi
cd ../

#Now f2py it!

# Intel Version
#f2py  --fcompiler=intel -c -m pyxlight_cs src_cs/pyxlight_cs.pyf src_cs/libxfoil_cs.a 

# gfortran Version
f2py  --fcompiler=gfortran -c -m pyxlight_cs src_cs/pyxlight_cs.pyf src_cs/libxfoil_cs.a 


#Test the Module
python module_import_test.py complex

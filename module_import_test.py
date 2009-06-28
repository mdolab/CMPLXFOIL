import sys,os

print 'Testing if module can be imported...'
print 

if len(sys.argv) > 1:
    if sys.argv[1] == 'real':
        real_test = True
        complex_test = False
    elif sys.argv[1] == 'complex':
        real_test = False
        complex_test = True
    else:
        real_test = True
        complex_test = True
else:
    real_test = True
    complex_test = True


if real_test:    
    try:
        import pyxlight
        print 'pyxlight.so imported successfully'
        os.system('mv pyxlight.so ./python')
        print 
    except:
        print 'Could not import pyxlight.so'
        print 
        sys.exit(1)
# end if

if complex_test:
    try:
        import pyxlight_cs
        print 'pyxlight_cs.so imported successfully'
        os.system('mv pyxlight_cs.so ./python')
        print 
    except:
        print 'Could not import pyxlight_cs.so'
        print 
        sys.exit(1)
# end if 

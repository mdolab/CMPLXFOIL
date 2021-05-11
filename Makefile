# Master makefile for pyXLIGHT. The actual makefile you want is:
# src/build/Makefile
# src_cs/build/Makefile

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. With the config file"; \
	echo " specified, rerun 'make' and the build will start"; \
	else make pyxlight_build;\
	fi;

clean:
	rm -fr src/build/*.mod
	rm -fr src/build/*.o
	rm -fr src/build/*.a
	rm -fr src/build/*.so
	rm -fr src_cs/build/*.mod
	rm -fr src_cs/build/*.o
	rm -fr src_cs/build/*.a
	rm -fr src_cs/build/*.so
	rm -f *~ config.mk;

pyxlight_build:
	ln -sf config/config.mk config.mk;
	(cd src/build/ && make)
	(cd src_cs/build/ && make)

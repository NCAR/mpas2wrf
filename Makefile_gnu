FC = gfortran
FFLAGS = -O2 -ffree-form -ffree-line-length-none
FCINCLUDES = $(shell nc-config --fflags)
FCLIBS = $(shell nc-config --flibs)

all:
	( cd src; $(MAKE) FC="$(FC)" FFLAGS="$(FFLAGS)" FCINCLUDES="$(FCINCLUDES)" FCLIBS="$(FCLIBS)" )
	if [ -e src/mpas2wrf_grid ] ; then \
	   ( cp src/mpas2wrf_grid . ) \
	fi;
clean:
	( cd src; $(MAKE) clean )
	rm -f mpas2wrf_grid 

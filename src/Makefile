.SUFFIXES: .F .o

OBJS = \
	scan_input.o \
	mpas_mesh.o \
	target_mesh.o \
	remapper.o \
	file_output.o \
	field_list.o \
	timer.o \
        constants_module.o \
        misc_definitions_module.o \
        module_map_utils.o \
        gridinfo_module.o \
        read_geogrid.o \
 	mpas2wrf_grid.o

all: $(OBJS)
	$(FC) -o mpas2wrf_grid $(OBJS) $(FCLIBS)

mpas2wrf_grid.o: scan_input.o mpas_mesh.o target_mesh.o file_output.o remapper.o field_list.o timer.o constants_module.o misc_definitions_module.o gridinfo_module.o module_map_utils.o read_geogrid.o 

remapper.o: mpas_mesh.o target_mesh.o scan_input.o

file_output.o: remapper.o scan_input.o

mpas_mesh.o: scan_input.o

field_list.o: scan_input.o

scan_input.o:

target_mesh.o:

gridinfo_module.o: constants_module.o misc_definitions_module.o scan_input.o 

module_map_utils.o: misc_definitions_module.o constants_module.o gridinfo_module.o remapper.o 



clean:
	rm -f *.mod *.o mpas2wrf_grid

.F.o:
	rm -f $@ $*.mod
	$(FC) $(FFLAGS) -c $*.F $(FCINCLUDES)

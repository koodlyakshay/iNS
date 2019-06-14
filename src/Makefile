PROG =	iNS

SRCS =	iNS.f90 mod_global.f90 \
	sub_compute_gradients.f90 sub_compute_massflux.f90 \
	sub_convective_residual.f90 \
	sub_flow_correction.f90 sub_linear_algebra.f90\
	sudb_output.f90 sub_preprocessing.f90 sub_pressure_residual.f90 \
	sub_read_input.f90 sub_time_integration.f90 sub_variable_operations.o sub_viscous_residual.f90 

OBJS =	iNS.o mod_global.o sub_compute_gradients.o \
	sub_compute_massflux.o sub_convective_residual.o \
	sub_flow_correction.o \
	sub_linear_algebra.o sub_output.o sub_preprocessing.o sub_pressure_residual.o \
	sub_read_input.o  sub_time_integration.o sub_variable_operations.o sub_viscous_residual.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

iNS.o: mod_global.o
sub_allocate_vars.o: mod_global.o
sub_compute_gradients.o: mod_global.o
sub_compute_massflux.o: mod_global.o
sub_convective_residual.o: mod_global.o
sub_flow_correction.o: mod_global.o
sub_initialize_vars.o: mod_global.o
sub_output.o: mod_global.o
sub_preprocessing.o: mod_global.o
sub_pressure_residual.o: mod_global.o
sub_read_input.o: mod_global.o
sub_time_integration.o: mod_global.o
sub_viscous_residual.o: mod_global.o

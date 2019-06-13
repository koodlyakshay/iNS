PROG =	iNS

SRCS =	iNS.f90 mod_global.f90 sub_allocate_vars.f90 \
	sub_compute_gradients.f90 sub_compute_massflux.f90 \
	sub_compute_pcgradients.f90 sub_convective_residual.f90 sub_explicit_euler.f90 \
	sub_gauss_seidel.f90 sub_implicit_euler.f90 sub_initialize_vars.f90 \
	sudb_output.f90 sub_preprocessing.f90 sub_pressure_residual.f90 \
	sub_read_input.f90 sub_viscous_residual.f90 

OBJS =	iNS.o mod_global.o sub_allocate_vars.o sub_compute_gradients.o \
	sub_compute_massflux.o sub_compute_pcgradients.o sub_convective_residual.o \
	sub_explicit_euler.o sub_gauss_seidel.o sub_implicit_euler.o \
	sub_initialize_vars.o sub_output.o sub_preprocessing.o sub_pressure_residual.o \
	sub_read_input.o  sub_viscous_residual.o

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
sub_compute_pcgradients.o: mod_global.o
sub_convective_residual.o: mod_global.o
sub_explicit_euler.o: mod_global.o
sub_implicit_euler.o: mod_global.o
sub_initialize_vars.o: mod_global.o
sub_output.o: mod_global.o
sub_preprocessing.o: mod_global.o
sub_pressure_residual.o: mod_global.o
sub_read_input.o: mod_global.o
sub_viscous_residual.o: mod_global.o

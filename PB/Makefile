PROG = iNS

SRCS =	iNS.f90 mod_global.f90 \
	sub_compute_gradients.f90

OBJS =	iNS.o iNS_backup.o iNS_upwind.o mod_global.o sub_compute_gradients.o

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
sub_compute_gradients.o: mod_global.o

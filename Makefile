FC=gfortran
F90FLAGS = -finit-local-zero -O1 -ffree-line-length-none -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

.f90.o:
	$(FC) $(F90FLAGS) -o $@ -c $<

OBJS = csfcal.o prep.o #detgen.o

all: csfcal

csfcal: $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

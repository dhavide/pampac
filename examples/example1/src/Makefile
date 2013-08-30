################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ../../../.
include ${TOPLEVEL}/Makefile.in  # Check Makefile.in for generic macros

F_SRCS = map_nu_omp.f coding.f
F90_SRCS = calculate_residual_f.f90 gmres.f90 single_corrector_step_f.f90
C_SRCS = C_interface.c main.c
OBJS = ${C_SRCS:.c=.o} ${F_SRCS:.f=.o} ${F90_SRCS:.f90=.o} 
LIBS += -llapack -lgfortran

################################################################################
# Targets
################################################################################

all: example1
	mv -i example1 ../

example1: $(OBJS)
	${MPICC} $(LIBPATHS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f ${OBJS} ../example1 *.d *~

################################################################################
# Dependencies and default compilation rules
################################################################################

# Awkward include/header dependencies for each object file
single_corrector_step_f.o: n_power.for constants.h
coding.o: n_power.for constants.h
map_nu_omp.o: set_parameters.f n_power.for parameters.h constants.h
calculate_residual_f.o: constants.h

################################################################################

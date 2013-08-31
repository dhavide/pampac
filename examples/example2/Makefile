################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ../..
include ${TOPLEVEL}/Makefile.in  # Parse Makefile.in for generic macros

NPROC ?= 4
OBJS = main.o compute_residual.o single_corrector_step.o setup_teardown.o
# Add dependence on libfftw3 into macro LIBS
# (macro LIBS is initialized in Makefile.in)
FFTW  = -lfftw3
LIBS += ${FFTW}
INCPATHS += -I /usr/include/atlas
################################################################################
# Targets
################################################################################
all: KSmain

exec: KSmain
	clear
	${MPIEXEC} -n $(NPROC) ./KSmain parameters.txt input.txt output.txt tree-data/tree

KSmain: $(OBJS)
	${MPICC} $(LIBPATHS) -o $@ $(OBJS) $(LIBS)

clean:
	/bin/rm -f *~ *.o *.d KSmain output.txt
	cd tree-data; make clean

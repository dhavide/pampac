################################################################################
# Macro definitions
################################################################################
MPICC   = /usr/bin/mpicc
MPIEXEC = /usr/bin/mpiexec
CFLAGS = -std=c99 -g -Wall

CSRC = $(wildcard *.c)
OBJS  = $(CSRC:.c=.o)


NPROC ?= 4
LIBPAMPAC=/home/dhavide/repositories/pampac/lib/libpampac.a
LIBPATHS = -L/usr/lib/atlas-base -L/home/dhavide/repositories/pampac/lib
INCPATHS = -I/home/dhavide/repositories/pampac/src

# Add dependence on libfftw3 into macro LIBS
# (macro LIBS is initialized in Makefile.in)

FFTW  = -lfftw3
MATH       = -lm
LAPACK     = /usr/lib/atlas-base/liblapack_atlas.so.3 /usr/lib/libcblas.so.3
BLAS       = -lgsl -lgslcblas
PAMPAC     = -lpampac
LIBS       = ${LAPACK} ${PAMPAC} ${FFTW}  ${BLAS} ${MATH}
INCPATHS += -I /usr/include/atlas
################################################################################
# Rules
################################################################################
.SUFFIXES: .o .c
.c.o:
	${MPICC} ${INCPATHS} ${CFLAGS} -c $< -o $@

################################################################################
# Targets
################################################################################
all: KSmain

exec: KSmain
	clear
	${MPIEXEC} -n $(NPROC) ./KSmain parameters.txt input.txt output.txt tree-data/tree

KSmain: $(OBJS)
	${MPICC} $(LIBPATHS) -o $@ $(OBJS) $(LIBS)

# /usr/lib/atlas-base/liblapack_atlas.so.3

clean:
	/bin/rm -f *~ *.o *.d KSmain output.txt
	cd tree-data; make clean

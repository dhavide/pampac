################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ..
#include ${TOPLEVEL}/Make.config  # Check Make.config for generic macros

MPICC = /usr/bin/gcc
CFLAGS = -std=c99 -g -I../src/ -I/usr/include/gsl/
#INCDIR = 
LDPATHS = -L../lib
LDLIBS = -lpampac -lgsl -lgslcblas -lm
OBJS  = 
NPROC ?=20
################################################################################
# Targets
################################################################################

svg:
	@cd data; make svg; cd ..

test_prune: test_setup.o test_prune.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)
test_construct: test_setup.o test_construct.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)

test_assign: test_setup.o test_assign.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)

test_choose: test_setup.o test_choose.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)

test_update: test_setup.o test_update.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)

test_repaint: test_setup.o test_repaint.o
	${MPICC} $(CFLAGS) $(LDPATHS)  -o $@ $^ $(LDLIBS)

assign: test_assign
	./$^ $(NPROC)

construct: test_construct
	./$^

choose: test_choose
	./$^

update: test_update
	./$^

repaint: test_repaint
	${MPIEXEC} -n 5 ./$^

prune: test_prune
	./$^

clean:
	rm -f test_construct test_update test_prune test_repaint test_assign test_choose *.o *.d data/*.gv data/*.svg *~ 

################################################################################

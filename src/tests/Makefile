################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ../..
include ${TOPLEVEL}/Makefile.in  # Check Makefile.in for generic macros

OBJS  = 
NPROC ?=20
################################################################################
# Targets
################################################################################
.PHONY: ../${LIBNAME}
	cd ..; make clean; make

test_prune: test_setup.o test_prune.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)
test_construct: test_setup.o test_construct.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)

test_assign: test_setup.o test_assign.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)

test_choose: test_setup.o test_choose.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)

test_update: test_setup.o test_update.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)

test_repaint: test_setup.o test_repaint.o
	${MPICC} $(LIBPATHS) -o $@ $^ $(LIBS)

execa: test_assign
	./$^ $(NPROC)
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 

exect: test_construct
	./$^
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 

execc: test_choose
	./$^
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 
	dot -Tsvg -o data/trees-00002.svg data/trees-00002.gv 
	dot -Tsvg -o data/trees-00003.svg data/trees-00003.gv 

execu: test_update
	./$^
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 

execr: test_repaint
	${MPIEXEC} -n 5 ./$^
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 

execp: test_prune
	./$^
	dot -Tsvg -o data/trees-00000.svg data/trees-00000.gv 
	dot -Tsvg -o data/trees-00001.svg data/trees-00001.gv 

clean:
	rm -f test_construct test_update test_prune test_repaint test_assign test_choose *.o *.d data/*.gv data/*.svg *~ 

################################################################################

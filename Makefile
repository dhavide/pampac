################################################################################
# Macro definitions
################################################################################

TOPLEVEL = .
include ${TOPLEVEL}/Makefile.in  # Parse Makefile.in for generic macros

################################################################################
# Targets
################################################################################
all:
	${MAKE} lib
	${MAKE} examples

lib:
	cd src; ${MAKE} lib

examples:
	cd examples; ${MAKE} all

clean:
	rm -f *~
	cd src; ${MAKE} clean
	cd examples; ${MAKE} clean

remake:
	${MAKE} clean
	${MAKE} all
################################################################################

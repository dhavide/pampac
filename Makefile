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
	${MAKE} example3

lib:
	cd src; ${MAKE} lib

example3:
	cd examples/example3; ${MAKE} KSmain

clean:
	rm -f *~
	cd src; ${MAKE} clean
	cd examples; ${MAKE} clean

remake:
	${MAKE} clean
	${MAKE} all   # NB: make all not done yet...
################################################################################

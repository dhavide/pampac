################################################################################
# Macro definitions
################################################################################
TOPLEVEL = ../.
include ${TOPLEVEL}/Makefile.in  # Check Makefile.in for generic macros
################################################################################
# Targets
################################################################################
all:
	${MAKE} example1
	${MAKE} example2

example1:
	cd example1; ${MAKE} all

example2:
	cd example2; ${MAKE} all
clean:
	cd example1; ${MAKE} clean
	cd example2; ${MAKE} clean


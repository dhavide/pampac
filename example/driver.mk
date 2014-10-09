TARGET := driver
SOURCES := main.c compute_residual.c helper_functions.c \
           single_corrector_step.c write_coordinates.c

# driver.config contains user specifications about the build & run
include ${EXAMPLEDIR}/driver.config

SUBMAKEFILES := ${EXAMPLEDIR}/trees.mk

# In TGT_POSTMAKE macro, redirection of stderr to /dev/null gets rid of an
# annoying error in moving the file driver on top of itself (which happens
# when "make driver" is invoked within the directory example).
TGT_POSTMAKE := mv -f driver ${EXAMPLEDIR}/driver 2> /dev/null; true
TGT_POSTCLEAN := rm -rf ./build

run: driver
	${MPIEXEC} -n ${NPROCS} ${EXAMPLEDIR}/driver ${PARAMETERFILE}

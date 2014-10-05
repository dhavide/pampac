TARGET := KS_example
SOURCES := \
   compute_residual.c helper_functions.c main.c single_corrector_step.c write_coordinates.c
TGT_PREREQS := lib/libpampac.a
TGT_POSTMAKE := mv KS_example example/KS_example
TGT_POSTCLEAN := rm -f example/KS_example
#SUBMAKEFILES = KS_exec.mk

# Modify these variables as appropriate for your system.
SRC_INCDIRS := ../src /usr/include/atlas
TGT_LDFLAGS := -Llib
TGT_LDLIBS := /usr/lib/atlas-base/liblapack_atlas.so.3 -lpampac -lgsl -lgslcblas -lm

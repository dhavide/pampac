TARGET := KS-example
SOURCES := \
   compute_residual.c \
   helper_functions.c \
   main.c \
   single_corrector_step.c
SRC_INCDIRS := ../src
TGT_PREREQS := lib/libpampac.a
SUBMAKEFILES = KS-exec.mk

# Modify these variables as appropriate for your system.
TGT_LDFLAGS := -Llib
TGT_LDLIBS := -llapacke -lpampac -lgsl -lgslcblas -lm

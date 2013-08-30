################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ..
include ${TOPLEVEL}/Makefile.in  # Check Makefile.in for generic macros

CSRC = $(wildcard *.c)
OBJS  = $(notdir $(CSRC:.c=.o))
# parse_options.o queue.o \
	  load_initial_coordinates.o  initialize_root_node.o compute_secant_direction.o \
	  init_PTnode.o delete_tree.o \
	  assign_processes.o prune_diverged_nodes.o count_children.o advance_root_node.o \
	  master_process.o slave_process.o  stop_slaves.o \
	  construct_predictor_nodes.o \
	  assess_residuals.o assign_predictor_steps.o compute_corrector_steps.o \
	  write_root_coordinates.o \
	  assign_color.o print_PTnode.o print_color.o print_tree.o visualize_tree.o \
	  construct_viable_paths.o choose_viable_paths.o

################################################################################
# Targets
################################################################################
all: ${OBJS}
	ar cvq ${LIBNAME} ${OBJS}

clean:
	rm -f ${OBJS} ${LIBNAME} *.d *~ 

################################################################################

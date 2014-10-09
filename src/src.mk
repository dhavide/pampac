TARGET := libpampac.a
SOURCES :=  advance_root_node.c assess_residuals.c \
            assign_depth.c assign_options.c \
            assign_predictor_steps.c assign_processes.c choose_viable_paths.c \
            compute_corrector_steps.c compute_secant_direction.c \
            construct_predictor_nodes.c construct_viable_paths.c \
            count_children.c create_root_node.c delete_tree.c \
            initialize_options.c initialize_secant.c init_PTnode.c \
            load_initial_coordinates.c master_process.c \
            parse_options.c principal_pampac_loop.c \
            print_state.c print_PTnode.c print_tree.c \
            prune_diverged_nodes.c queue.c slave_process.c \
            stop_slaves.c validate_options.c visualize_tree.c write_root_coordinates.c

TGT_POSTMAKE := mkdir -p ${TOPDIR}/lib; mv libpampac.a ${TOPDIR}/lib
TGT_POSTCLEAN := rm -f ${TOPDIR}/lib/libpampac.a; rm -rf ${TOPDIR}/src/build/


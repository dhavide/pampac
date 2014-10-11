#include "mpi.h"
#include "pampac.h"

void
principal_pampac_loop (int N_p, PTnode * root, options_struct * opts) {
  double time_init, time_final, lambda, lambda_min, lambda_max, h, h_min;
  int lambda_index, global_iter, max_global_iter;
  bool has_failed, has_completed;

  if (opts->verbose>3)
    printf ("principal_pampac_loop: Beginning main computation...\n");

  /* Setting convenient aliases for optional parameters */
  lambda_min = opts->lambda_min;
  lambda_max = opts->lambda_max;
  lambda_index = opts->lambda_index;
  lambda = root->z[lambda_index];
  h_min = opts->h_min;
  h = root->h;
  max_global_iter = opts->max_global_iter;
  if (opts->verbose > 0) {
    printf ("principal_pampac_loop: h=%7.1e, h_min=%7.1e\n", h, h_min);
    printf ("principal_pampac_loop: Maximum global iterations=%d\n",
            max_global_iter);
  }

  global_iter = 0;
  has_failed = false;
  has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
  if (opts->verbose > 0)
    printf ("principal_pampac_loop: Beginning global iteration.\n");

  /* Log image of initial node if necessary */
  if (opts->verbose > 2)
    visualize_tree (root, opts);

  time_init = MPI_Wtime ();
  while (!has_completed && !has_failed) {
    printf ("principal_pampac_loop: while iteration %i:\n", global_iter);
    has_failed = (global_iter > max_global_iter) || (h < h_min);
    if (has_failed) {
      printf ("principal_pampac_loop: Premature termination\n"); 
      printf ("(global_iter = %i > %i = max_global_iter) or ",
               global_iter, max_global_iter);
      printf ("(h = %7.1e < %7.1e = h_min)\n", h, h_min);
      break;
    }
    /* Spawn new nodes on tree at leaves if possible. */
    construct_predictor_nodes (root, opts);
    assign_processes (root, N_p);
    assign_predictor_steps (root, opts);
    visualize_tree (root, opts);

    compute_corrector_steps (root, N_p);
    assess_residuals (root, opts);
    visualize_tree (root, opts);

    prune_diverged_nodes (root, opts);
    visualize_tree (root, opts);

    construct_viable_paths (root);
    choose_viable_paths (root);
    visualize_tree (root, opts);

    advance_root_node (&root, opts);
    visualize_tree (root, opts);

    global_iter++;
    lambda = root->z[lambda_index];
    has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
    h = root->h;
  }
  time_final = MPI_Wtime ();

  if (opts->verbose > 0) {
    if (has_completed) {
      printf ("principal_pampac_loop: Completed main loop\n");
      printf ("principal_pampac_loop: %7.1e <= lambda = %7.1e <= %7.1e\n",
              lambda_min, lambda, lambda_max); 
    }
    printf ("principal_pampac_loop: Time elapsed = %g\n", time_final - time_init);
    printf ("principal_pampac_loop: global iterations = %d.\n", global_iter);
  }
  return;
}

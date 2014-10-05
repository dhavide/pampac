#include "mpi.h"
#include "pampac.h"

void
principal_pampac_loop (int N_p, PTnode * root, options_struct * opts) {
  double time_init, time_final, lambda, lambda_min, lambda_max, h, h_min;
  int lambda_index, global_iter, max_global_iter;
  bool has_failed, has_completed;

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
    has_failed = (global_iter > max_global_iter) || (h < h_min);
    if (has_failed) {
      printf ("principal_pampac_loop:");
      printf (" Failed to attain goal\n");
      printf ("global iterations=%d, h=%10.5g\n", global_iter, h);
      break;
    }
    /* Spawn new nodes on tree at leaves if possible. */
    construct_predictor_nodes (root, opts);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Constructed predictor nodes...\n");
    assign_processes (root, N_p);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Assigned processes...\n");
    assign_predictor_steps (root, opts);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Assigned predictor steps...\n");
    if (opts->verbose > 2)
      visualize_tree (root, opts);

    compute_corrector_steps (root, N_p);

    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Computing (concurrent) corrector steps...\n");

    assess_residuals (root, opts);

    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Assessed residuals...\n");

    if (opts->verbose > 2)
      visualize_tree (root, opts);

    prune_diverged_nodes (root, opts);

    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Pruned diverged nodes...\n");

    if (opts->verbose > 2)
      visualize_tree (root, opts);

    construct_viable_paths (root);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Constructed viable paths...\n");

    choose_viable_paths (root);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Choose viable paths...\n");

    if (opts->verbose > 2)
      visualize_tree (root, opts);

    advance_root_node (&root, opts);
    if (opts->verbose > 3)
      printf ("principal_pampac_loop: Advanced root node...\n");

    if (opts->verbose > 2)
      visualize_tree (root, opts);

    global_iter++;
    lambda = root->z[lambda_index];
    has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
    h = root->h;
  }
  time_final = MPI_Wtime ();

  if (opts->verbose > 0) {
    if (has_completed)
      printf ("principal_pampac_loop: Completed main loop\n");
    else
      printf ("principal_pampac_loop: Premature termination\n");
    printf ("Time elapsed: %g\n", time_final - time_init);
    printf ("global iterations=%d.\n", global_iter);
    printf ("has_completed = %5s, has_failed = %5s\n",
            (has_completed ? " true" : "false"),
            (has_failed ? " true" : "false"));
    printf ("lambda =%10g, lambda_min=%10g, lambda_max=%10g\n",
            lambda, lambda_min, lambda_max);
    printf ("h = %10g, h_min = %10g\n", h, h_min);
  }
  return;
}

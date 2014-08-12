#include "mpi.h"
#include "pampac.h"

void principal_pampac_loop (int N_p, PTnode *root, options_struct *opts)
{
  double time_init, time_final, lambda, lambda_min, lambda_max,
         h, h_min, h_max;
  int lambda_index, global_iter, max_global_iter;
  bool has_failed, has_completed;

  /* Setting convenient aliases for optional parameters */
  lambda_min = opts->lambda_min;
  lambda_max = opts->lambda_max;
  lambda_index = opts->lambda_index;
  lambda = root->z[lambda_index];
  h_min = opts->h_min;
  h_max = opts->h_max;
  max_global_iter = opts->max_global_iter;

  //tree_file_count = 0;
  global_iter = 0;
  has_failed = false;
  has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
  if (opts->verbose>0)
    printf ("principal_pampac_loop: Beginning global iteration.\n");

//               visualize_tree (root, opts, tree_file_count++)
//  LOG_TREE(NULL,1); /* Log image of initial node if necessary */
  time_init = MPI_Wtime ();
  while (!has_completed && !has_failed)
    {
      has_failed = (global_iter>max_global_iter) || (h<h_min);
      if (has_failed)
      {
          printf("principal_pampac_loop:");
          printf(" Failed to attain goal\n");
          break;
      }
      /* Spawn new nodes on tree at leaves if possible. */
      construct_predictor_nodes (root, opts);
      assign_processes (root, N_p);
assign_predictor_steps (root, opts);
      compute_corrector_steps (root, N_p);
assess_residuals (root, opts);
prune_diverged_nodes (root, opts);
      construct_viable_paths (root);
choose_viable_paths (root);
advance_root_node (&root, opts);
      global_iter++;
      lambda = root->z[lambda_index];
      has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
      h = root->h;
    }
  time_final = MPI_Wtime ();

  if (opts->verbose>0)
    {
      if (has_completed)
    {
      printf ("principal_pampac_loop: completed main loop\n");
      printf ("Time elapsed: %g\n", time_final - time_init);
      printf ("Global_Iter = %d global iterations.\n", global_iter);
      printf ("has_completed = %5s, has_failed = %5s\n",
              (has_completed ? " true" : "false"),
              (has_failed ? " true" : "false"));
      printf ("lambda =%10g, lambda_min=%10g, lambda_max=%10g\n",
             lambda, lambda_min, lambda_max);
      printf ("h = %10g, h_min = %10g, h_max = %10g\n",h,h_min,h_max);
    }
      else
      {
    printf ("Premature termination: ");
    printf ("has_completed = %5s, has_failed = %5s\n",
              (has_completed ? " true" : "false"),
              (has_failed ? " true" : "false"));
    }
    }
return;
}

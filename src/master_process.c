#include "mpi.h"
#include "pampac.h"
/**********************************************************************/
/* Main routine executed by master processor. The bulk of the work is */
/* in the routine construct_continuation_curve; the remaining code is */
/* largely for preparation and clean-up.                              */
/**********************************************************************/
void
master_process (int p_id, int N_p, options_struct * opts)
{
  double time_init, time_final, lambda, lambda_min, lambda_max, 
         h, h_min, h_max;
  int lambda_index, tree_file_count, global_iter, max_global_iter;
  bool has_failed, has_completed;
  PTnode *root = NULL;
  initialize_root_node (&root, opts);
  lambda_min = opts->lambda_min;
  lambda_max = opts->lambda_max;
  lambda_index = opts->lambda_index;
  h_min = opts->h_min;
  h_max = opts->h_max;
  lambda = root->z[lambda_index];
  max_global_iter = opts->max_global_iter;

  tree_file_count = 0;
  global_iter = 0;
  has_failed = false;
  has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
  if (opts->verbose>0)
    printf ("master_process: beginning global iteration.\n");
  time_init = MPI_Wtime ();

  /* Macro to execute command and log image of tree if requested */
#define LOG_TREE(comm,n) (comm);		  \
                     if (opts->verbose>n) \
		       visualize_tree (root, opts, tree_file_count++)
  LOG_TREE(NULL,1); /* Log image of initial node if necessary */
  while (!has_completed && !has_failed)
    {
      /* Spawn new nodes on tree at leaves if possible. */
      construct_predictor_nodes (root, opts);
      assign_processes (root, N_p);
      LOG_TREE(assign_predictor_steps (root, opts),1);
      compute_corrector_steps (root, N_p);
      LOG_TREE(assess_residuals (root, opts),1);
      LOG_TREE(prune_diverged_nodes (root, opts),1);
      construct_viable_paths (root);
      LOG_TREE(choose_viable_paths (root),1);
      LOG_TREE(advance_root_node (&root, opts),1);
      global_iter++;
      lambda = root->z[lambda_index];
      has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
      h = root->h;
      has_failed = (global_iter>max_global_iter) || (h<h_min);
    }
  time_final = MPI_Wtime ();
  if  (opts->verbose>0)
    {
      if (has_completed)
	{
	  printf ("master_process: completed main loop\n");
	  printf ("Time elapsed: %g\n", time_final - time_init);
	  printf ("Global_Iter = %d global iterations.\n", global_iter);
	  printf ("has_completed = %d, has_failed = %d\n", has_completed, has_failed);
	  printf ("lambda =%10g, lambda_min=%10g, lambda_max=%10g\n",
		  lambda, lambda_min, lambda_max);
	  printf ("h = %10g, h_min = %10g, h_max = %10g\n",h,h_min,h_max);
	}
      else
	printf ("Premature termination:has_completed = %d, has_failed = %d\n",
		has_completed, has_failed);
    }
  write_root_coordinates (root, opts);
  stop_slaves(N_p);
  free (root->z_parent);
  root->z_parent = NULL;
  delete_tree (root);
  root = NULL;
  return;
}

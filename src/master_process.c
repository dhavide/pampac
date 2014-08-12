#include "mpi.h"
#include "pampac.h"
#include <gsl/gsl_cblas.h>
/**********************************************************************/
/* Main routine executed by master processor. The bulk of the work is */
/* in the routine principal_pampac_loop; the remaining code is        */
/* largely for preparation and clean-up.                              */
/**********************************************************************/
void
master_process (int N_p, options_struct *opts) {
  PTnode *root = NULL;
  double r_nrm;
  int lambda_index = opts->lambda_index;
  bool has_succeeded;

  /* First verification: meaningful tree depth must be larger than 1 */
  if (opts->max_depth<=1) {
    printf ("master_process: Require max_depth > 1; max_depth = %d\n",
            opts->max_depth);
    printf ("Terminating...\n");
    goto cleanup;
  }

  /* Second verification: Allocate root_node successfully */
  has_succeeded = create_root_node (&root, opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to allocate memory for root node.\n");
    printf ("Terminating...\n");
    goto cleanup;
  }

  /* Third verification: Load initial point from disk */
  has_succeeded = load_initial_coordinates (root, opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to read first point into root node.\n");
    printf ("Terminating...\n");
    goto cleanup;
  }

  if (opts->verbose>0) {
    printf ("master_process: Loaded first point from file.\n");
    compute_residual (root->N_dim, root->z, root->T_init);
    r_nrm = cblas_dnrm2 (root->N_dim, root->T_init, 1);
    printf ("Initial residual = %12.5g\n", r_nrm);
    printf ("Initial lambda = %12.5g\n", root->z[lambda_index]);
  }

  /* Fourth verification: Determine initial secant direction */
  has_succeeded = initialize_secant (root, opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to determine initial secant direction.\n");
    printf ("Maximum of %d iterations exceeded.\n", opts->max_iter);
    printf ("Desired residual tolerance: %12.5e\n",opts->tol_residual);
    printf ("Terminating.\n");
    goto cleanup;
  }

  principal_pampac_loop (N_p, root, opts);
  write_root_coordinates (root, opts);

cleanup:
  if (opts->verbose>0)
    printf("master_process: Shutting down slave processes.\n");
  stop_slaves(N_p);
  if (opts->verbose>0)
    printf("master_process: Cleaning up memory.\n");
  free (root->z_init);
  root->z_init = NULL;
  delete_tree (root);
  root = NULL;
  return;
}

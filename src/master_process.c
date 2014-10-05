#include "mpi.h"
#include "pampac.h"
#include <gsl_cblas.h>
/**********************************************************************/
/* Main routine executed by master processor. The bulk of the work is */
/* in the routine principal_pampac_loop; the remaining code is        */
/* largely for preparation and clean-up.                              */
/**********************************************************************/
void
master_process (int N_p, int argc, char *argv[]) {
  PTnode *root = NULL;
  options_struct opts;
  bool has_succeeded;

  /* Default options before passing */
  initialize_options (&opts);
  opts.verbose = 2;
  opts.max_children = 0;
  has_succeeded = parse_options (argc, argv, &opts);
  /* Verification: meaningful tree depth must be larger than 1 */
  if (!has_succeeded) {
    printf ("master_process: Unable to parse parameter file\n");
    printf ("Terminating...\n");
    goto cleanup;
  }

  /* Verification: meaningful tree depth must be larger than 1 */
  if (opts.max_depth<=1) {
    printf ("master_process: Require max_depth > 1; max_depth = %d\n",
            opts.max_depth);
    printf ("Terminating...\n");
    goto cleanup;
  }

  /* Verification: Allocate root_node successfully */
  has_succeeded = create_root_node (&root, &opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to allocate memory for root node.\n");
    printf ("Terminating...\n");
    goto cleanup;
  }

  /* Verification: Load initial point from disk */
  has_succeeded = load_initial_coordinates (root, &opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to read first point into root node.\n");
    printf ("Terminating...\n");
    goto cleanup;
  }

  if (opts.verbose>0) {
    double res_nrm;
    int lambda_index = opts.lambda_index;
    printf ("master_process: Loaded first point from file.\n");
    compute_residual (root->N_dim, root->z, root->T_init);
    res_nrm = cblas_dnrm2 (root->N_dim-1, root->T_init, 1);
    printf ("master_process: Initial residual = %12.5g\n", res_nrm);
    printf ("master_process: ");
    printf ("Initial lambda = %12.5g\n", root->z[lambda_index]);
  }

  /* Verification: Determine initial secant direction */
  has_succeeded = initialize_secant (root, &opts);
  if (!has_succeeded) {
    printf ("master_process:");
    printf (" Failed to determine initial secant direction.\n");
    printf ("Maximum of %d iterations exceeded.\n", opts.max_iter);
    printf ("Desired residual tolerance: %12.5e\n",opts.tol_residual);
    printf ("Terminating.\n");
    goto cleanup;
  }

  principal_pampac_loop (N_p, root, &opts);
  write_root_coordinates (root, &opts);

cleanup:
  if (opts.verbose>0)
    printf("master_process: Shutting down slave processes.\n");
  stop_slaves(N_p);
  if (opts.verbose>0)
    printf("master_process: Cleaning up memory.\n");
  if (opts.input_filename!=NULL)
    free (opts.input_filename);
  if (opts.output_filename!=NULL)
    free (opts.output_filename);
  if (opts.scale_factors!=NULL)
    free (opts.scale_factors);
  if (root!=NULL) {
    /* z_init is a shallow copy in every node of the tree *except* the
    * root node. As such, z_init must be *explicitly* deallocated
    * prior to calling the function delete_tree. */
    if (root->z_init != NULL)
      free (root->z_init);
    delete_tree (root);
  }
  return;
}

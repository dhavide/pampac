#include <gsl/gsl_cblas.h>
#include "pampac.h"
void initialize_root_node (PTnode** root_addr, options_struct *opts)
{
  PTnode *root;
  int k, N_dim = opts->N_dim;
  if (opts->verbose>0)
    printf("initialize_root_node: setting fields of root node.\n");
  /* Fill in fields of root node z appropriately */
  root = init_PTnode (opts->max_children);
  root->N_dim = N_dim;
  root->depth = 0;
  root->color = GREEN;
  root->nu = 0;
  root->h = opts->h_init;
  root->z = malloc (N_dim * sizeof (double));
  root->z_init = malloc (N_dim * sizeof (double));
  root->T_init = malloc (N_dim * sizeof (double));
  load_initial_coordinates (root, opts);

  if (opts->verbose>0)
  {
    printf ("initialize_root_node: Loaded first point from file.\n");
    compute_residual (root->N_dim, root->z, root->T_init);
    print_PTnode (root);
    double r_nrm = cblas_dnrm2 (root->N_dim, root->T_init, 1);
    printf("Initial residual = %g\n", r_nrm);
    printf("Initial lambda = %g\n", root->z[opts->lambda_index]);
  }

  for (k=0; k<N_dim; k++)
    root->z_init[k] = root->z[k];

  root->z_init[opts->lambda_index] -= (opts->lambda_dir)*(opts->delta_lambda);
  /* To find another point on the continuation curve, search along the */
  /* drection of a unit vector with all zeros except (-1.0) in the     */
  /* component of lambda_index (index of continuation parameter).      */
  for (k = 0; k < N_dim - 1; k++)
    root->T_init[k] = 0.e0;
  /* Secant is elementary unit vector in reverse direction. */
  root->T_init[opts->lambda_index] = 1.0;

  bool has_converged=false, has_failed=false;
  int count = 0;
  double residual[N_dim-1], r_nrm;
  while ((!has_converged) && (!has_failed))
    {
      count++;
      single_corrector_step (N_dim, root->z_init, root->T_init);
      compute_residual (N_dim, root->z_init, residual);
      r_nrm = cblas_dnrm2 (N_dim, residual, 1);
      if (opts->verbose>1)
          printf("initialize_root_node: count=%d, r_nrm=%g.\n",
                  count, r_nrm);
      has_converged = (r_nrm < opts->tol_residual);
      has_failed = (count > opts->max_iter);
    }
  if (has_failed)
    {
      printf("Failed to determine second point on continuation curve.\n");
      printf("Maximum of %d corrector iterations exceeded.\n", opts->max_iter);
      printf("Desired residual tolerance: %12.5e\n",opts->tol_residual);
      printf("Aborting processes.\n");
      return;
    }
  compute_secant_direction (root);
  write_root_coordinates (root, opts);
  *root_addr = root;
  return;
}

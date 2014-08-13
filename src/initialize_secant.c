#include <gsl/gsl_cblas.h>
#include "pampac.h"
bool initialize_secant (PTnode* root, options_struct *opts) {
  int count, k, N_dim;
  bool has_converged, has_failed;
  double *residual, r_nrm;

  /* Assume that node root has been initialized with a vector
   * root->z somehow. This function determines a nearby point on
   * the homotopy curve by perturbing the continuation parameter
   * a small amount and using corrector steps to return to the
   * curve. The two nearby points determine a secant direction. */
  N_dim = opts->N_dim;
  residual = malloc((N_dim-1)*sizeof(double));
  /* Copy z into z_init */
  for (k=0; k<N_dim; k++)
    root->z_init[k] = root->z[k];

  /* Modify z_init by small perturbation in direction lambda_index */
  root->z_init[opts->lambda_index] -= (opts->lambda_dir) *
                                      (opts->delta_lambda);
  /* Initial secant direction is coordinate vector in the direction
   * of lambda_index (index of continuation parameter).  */
  for (k = 0; k < N_dim - 1; k++)
    root->T_init[k] = 0.e0;
  root->T_init[opts->lambda_index] = 1.0;

  if (opts->verbose>0) {
    printf("initialize_secant: Iteration to get initial");
    printf(" secant direction.\n");
  }
  count = 0;
  has_converged = false;
  has_failed = false;
  while (!has_converged) {
    count++;
    has_failed = (count > opts->max_iter);
    if (has_failed) {
      printf("initialize_secant: Failed to determine initial secant");
      printf(" direction.\n");
      printf("Maximum of %d corrector iterations attained.\n",
             opts->max_iter);
      printf("Residual norm: %12.5e\n",opts->tol_residual);
      printf("Desired residual norm: %12.5e\n",opts->tol_residual);
      printf("Aborting processes.\n");
      return false;
    }
    single_corrector_step (N_dim, root->z_init, root->T_init);
    compute_residual (N_dim, root->z_init, residual);
    r_nrm = cblas_dnrm2 (N_dim-1, residual, 1);
    if (opts->verbose>1) {
      printf("initialize_secant: count=%3d,", count);
      printf(" residual norm=%12.5g.\n", r_nrm);
    }
    has_converged = (r_nrm < opts->tol_residual);
  }
  compute_secant_direction (root);
  return true;
}

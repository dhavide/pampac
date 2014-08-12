#include <mpi.h>
#include "pampac.h"
/**********************************************************************/
/* This function traverses the tree, and, for each node in the tree,  */
/* receives the residual as computed by one of the processors. The    */
/* value received is used to assign colors to the nodes which is the  */
/* basis of the decision-making algorithm for advancing on the curve. */
/**********************************************************************/
void
assess_residuals (PTnode * node, options_struct * opts) {
  int k, n_received;
  double residual_old, TOL, GAMMA, MU;
  bool has_failed, has_converged, has_almost_converged;
  MPI_Status status;
  TOL = opts->tol_residual;
  GAMMA = opts->gamma;
  MU = opts->mu;

  for (k = 0; k < node->max_children; k++)
    if (node->child[k] != NULL)
      assess_residuals (node->child[k], opts); /* Recursive call */

  if ((node->color==GREEN) || (node->pid<0))
    return;

  residual_old = node->res_norm; /* For measuring progress later */
  /* Corresponding MPI_Send calls in slave_process.
     Note: if node->pid==MPI_PROC_NULL, MPI_Recv will not hang. */

  MPI_Recv (node->z, node->N_dim, MPI_DOUBLE, node->pid,
            CONTINUE_TAG, MPI_COMM_WORLD, &status);
  /* Verifies that data was actually received */
  MPI_Get_count (&status, MPI_DOUBLE, &n_received);
  if (n_received==node->N_dim)
    node->nu++;
  MPI_Recv (&(node->res_norm), 1, MPI_DOUBLE, node->pid,
            CONTINUE_TAG, MPI_COMM_WORLD, &status);

  /* RED    -> keep computing corrector steps (default)
     GREEN  -> converged; no more corrector steps needed
     YELLOW -> almost converged; one more corrector step
     BLACK  -> insufficient progress or diverged         */

  has_failed = (node->nu > opts->max_iter) ||
               (log10(node->res_norm) > log10(MU) + log10(residual_old));
  has_converged = (node->res_norm < TOL);
  has_almost_converged = (GAMMA * log10 (node->res_norm) < log10 (TOL));
  node->color = RED;
  if (has_converged)
    node->color = GREEN;
  else if (has_failed)
    node->color = BLACK;
  else if (has_almost_converged)
    node->color = YELLOW;
  if (opts->verbose>0) {
    printf ("Pre/Post-residual=%10.4e/%-10.4e, nu=%d, h=%5g, color=",
            residual_old, node->res_norm, node->nu, node->h);
    print_color(node,stdout);
    printf("\n");
  }
  return;
}


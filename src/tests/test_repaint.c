#include "mpi.h"
#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main (int argc, char *argv[])
{
  int p_id, N_p;
  MPI_Init (&argc, &argv);	           /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id);   /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);	   /* find total number of processes */

  if (p_id==0) {
  /* Constructs wide string of nodes to test assess_residuals. */
  options_struct opts = set_options (/*depth=*/2, /*width=*/5);
  opts.tol_residual = 1.0e-6;
  PTnode* root = make_tree (&opts);

  /* Assign values to node fields. */
  for (int k=0; k<opts.max_children; k++)
    {
      root->child[k]->label = k+1;
      root->child[k]->pid = -1;
    }
  root->res_norm = 1.0e-7;
  root->child[4]->pid = 1; // To actually receive message
  root->child[0]->res_norm = 0.99e-4;
  root->child[1]->res_norm = 1.01e-4;
  root->child[2]->res_norm = 0.99e-6;
  root->child[3]->res_norm = 0.5e-4;
  root->child[4]->res_norm = 1.0e-4;

  int count=0;
  visualize_tree(root,&opts,count++);
  assess_residuals (root, &opts);
  print_tree (root);
  visualize_tree(root,&opts,count++);
  delete_tree(root);
  } else if (p_id==1) {
    double my_vals[6] = {0.1,-0.2,0.3,-0.4,1.0e2,-1.1e1};
    double r_norm = 1.45e-4;
    MPI_Send (my_vals, 6, MPI_DOUBLE, 0, CONTINUE_TAG, MPI_COMM_WORLD);
    MPI_Send (&r_norm, 1, MPI_DOUBLE, 0, CONTINUE_TAG, MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

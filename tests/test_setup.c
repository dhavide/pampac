#include "pampac.h"

/* A sample tree to use for testing pampac library routines. */

options_struct
set_options(int depth, int max_children) {
  options_struct opts;
  int k;
  opts.N_dim = 6;
  opts.lambda_min = 0.0;
  opts.lambda_max = 1.0;
  opts.lambda_index = opts.N_dim-1;
  opts.lambda_dir = 1;
  opts.delta_lambda = 1.0e-5;
  opts.h_min = 1.0e-10;
  opts.h_max = 1.0e2;
  opts.h_init = 1.0e-1;
  opts.max_iter = 5;
  opts.tol_residual = 1.0e-6;
  opts.mu = 1.3;
  opts.gamma = 1.5;
  opts.max_depth = depth;
  opts.max_children = max_children;
  opts.scale_process[0] = 0.5;
  for (k=0; k<max_children; k++)
    opts.scale_process[k] = 0.5 + k*0.25;
  opts.output_file_name = "foo.txt";
  opts.input_file_name = "foo.txt";
  opts.base_name_tree = "./data/trees";
  return (opts);
}

PTnode* make_tree(options_struct *opts) {
  int k;
  PTnode *root = init_PTnode (opts->max_children);
  root->N_dim = opts->N_dim;
  root->color = GREEN;
  root->h = opts->h_init;
  root->nu = 0;
  root->max_children = opts->max_children;
  root->depth = 0;
  root->z = malloc (root->N_dim * sizeof (double));
  root->T_parent = malloc (root->N_dim * sizeof (double));
  root->z_parent = malloc (root->N_dim * sizeof (double));
  for (k=0; k<root->N_dim; k++) {
    root->z[k] = (double)(k+2);
    root->z_parent[k] = (double)(k+1);
  }
  compute_secant_direction (root);
  /* Generate new levels of tree (nodes only) */
  for (k=1; k<opts->max_depth; k++) {
    construct_predictor_nodes (root, opts);
    assign_processes (root, 50);  /* Assume large number of processors */
    assign_predictor_steps (root, opts);
  }
  return (root);
}


void tree43(PTnode* root) {
  /* Given tree of depth 4 with 3 children per node, assigns various
     colors to the nodes and prunes a few subtrees. */
  root->color = GREEN;
  root->nu = 4;

  root->child[0]->color = GREEN;
  root->child[1]->color = GREEN;
  root->child[2]->color = GREEN;

  root->child[0]->nu = 4;
  root->child[1]->nu = 5;
  root->child[2]->nu = 3;

  root->child[0]->child[0]->color = GREEN;
  root->child[0]->child[1]->color = BLACK;
  root->child[0]->child[2]->color = YELLOW;

  root->child[0]->child[0]->nu = 4;
  root->child[0]->child[1]->nu = 3;
  root->child[0]->child[2]->nu = 3;

  root->child[1]->child[0]->color = BLACK;
  delete_tree (root->child[1]->child[1]);
  root->child[1]->child[1]=NULL;
  root->child[1]->child[2]->color = BLACK;

  root->child[1]->child[0]->nu = 4;
  root->child[1]->child[2]->nu = 3;

  root->child[2]->child[0]->color = YELLOW;
  root->child[2]->child[1]->color = RED;
  root->child[2]->child[2]->color = BLACK;

  root->child[2]->child[0]->nu = 4;
  root->child[2]->child[1]->nu = 2;
  root->child[2]->child[2]->nu = 3;

  root->child[0]->child[0]->child[0]->color = BLACK;
  delete_tree (root->child[0]->child[0]->child[1]);
  root->child[0]->child[0]->child[1]=NULL;
  root->child[0]->child[0]->child[2]->color = BLACK;

  root->child[0]->child[0]->child[0]->nu = 3;
  root->child[0]->child[0]->child[2]->nu = 3;

  root->child[0]->child[1]->child[0]->color = RED;
  root->child[0]->child[1]->child[1]->color = YELLOW;
  root->child[0]->child[1]->child[2]->color = RED;

  root->child[0]->child[1]->child[0]->nu = 3;
  root->child[0]->child[1]->child[1]->nu = 3;
  root->child[0]->child[1]->child[2]->nu = 3;

  root->child[0]->child[2]->child[0]->color = RED;
  root->child[0]->child[2]->child[1]->color = BLACK;
  root->child[0]->child[2]->child[2]->color = YELLOW;

  root->child[0]->child[2]->child[0]->nu = 2;
  root->child[0]->child[2]->child[1]->nu = 2;
  root->child[0]->child[2]->child[2]->nu = 2;

  root->child[1]->child[0]->child[0]->color = RED;
  root->child[1]->child[0]->child[1]->color = YELLOW;
  root->child[1]->child[0]->child[2]->color = RED;

  root->child[1]->child[0]->child[0]->nu = 2;
  root->child[1]->child[0]->child[1]->nu = 4;
  root->child[1]->child[0]->child[2]->nu = 2;

  root->child[1]->child[2]->child[0]->color = RED;
  root->child[1]->child[2]->child[1]->color = GREEN;
  root->child[1]->child[2]->child[2]->color = YELLOW;

  root->child[1]->child[2]->child[0]->nu = 3;
  root->child[1]->child[2]->child[1]->nu = 5;
  root->child[1]->child[2]->child[2]->nu = 2;

  root->child[2]->child[0]->child[0]->color = GREEN;
  root->child[2]->child[0]->child[1]->color = YELLOW;
  root->child[2]->child[0]->child[2]->color = RED;

  root->child[2]->child[0]->child[0]->nu = 5;
  root->child[2]->child[0]->child[1]->nu = 2;
  root->child[2]->child[0]->child[2]->nu = 2;

  root->child[2]->child[1]->child[0]->color = BLACK;
  root->child[2]->child[1]->child[1]->color = GREEN;
  delete_tree (root->child[2]->child[1]->child[2]);
  root->child[2]->child[1]->child[2] = NULL;

  root->child[0]->child[2]->child[0]->nu = 2;
  root->child[0]->child[2]->child[1]->nu = 5;

  root->child[2]->child[2]->child[0]->color = YELLOW;
  root->child[2]->child[2]->child[1]->color = GREEN;
  root->child[2]->child[2]->child[2]->color = BLACK;

  root->child[2]->child[2]->child[0]->nu = 4;
  root->child[2]->child[2]->child[0]->nu = 4;
  root->child[2]->child[2]->child[0]->nu = 2;

  return;
}

#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main (int argc, char* argv[]) {
  options_struct opts = set_options (1,2);
  PTnode* root = make_tree (&opts);
  int nproc = atoi(argv[1]), count=0;
  opts.max_depth = 4;
  construct_predictor_nodes (root,&opts);
  assign_processes(root, nproc);
  assign_predictor_steps (root, &opts);
  visualize_tree(root,&opts, count++);
  print_tree (root);
  printf("*********************************************************\n");
  construct_predictor_nodes (root,&opts);
  assign_processes(root, nproc);
  assign_predictor_steps (root, &opts);
  delete_tree (root->child[0]->child[1]);
  root->child[0]->child[1] = NULL;
  delete_tree (root->child[1]->child[0]);
  root->child[1]->child[0] = NULL;
  delete_tree (root->child[1]->child[1]);
  root->child[1]->child[1] = NULL;
  printf("*********************************************************\n");
  construct_predictor_nodes (root,&opts);
  assign_processes(root, nproc);
  assign_predictor_steps (root, &opts);
  visualize_tree(root,&opts, count++);
  print_tree (root);
  printf("Calling delete_tree...\n");
  free (root->z_parent); /* Only one needs deallocation */
  delete_tree(root);
  return 0;
}

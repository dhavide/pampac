#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main ()
{
  options_struct opts = set_options (1,2);
  PTnode* root = make_tree (&opts);
  opts.max_depth = 3;
  construct_predictor_nodes (root,&opts);
  assign_processes(root, 20);
  construct_predictor_nodes (root,&opts);
  assign_processes(root, 20);
  delete_tree (root->child[0]->child[1]);
  root->child[0]->child[1] = NULL;
  delete_tree (root->child[1]->child[0]);
  root->child[1]->child[0] = NULL;
  delete_tree (root->child[1]->child[1]);
  root->child[1]->child[1] = NULL;
  int count=0;
  visualize_tree(root,&opts,count++);
  construct_predictor_nodes (root,&opts);
  assign_processes(root, 20);
  visualize_tree(root,&opts,count++);
  printf("Calling delete_tree...\n");
  free (root->z_parent); /* Only one needs deallocation */
  delete_tree(root);
  return 0;
}

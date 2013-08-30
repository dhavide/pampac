#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main ()
{
  options_struct opts = set_options (3,2);
  PTnode* root = make_tree (&opts);
  /* Assign values to node fields. */
  root->nu = 2;
  root->child[0]->color = GREEN;
  root->child[0]-> nu = 3;
  root->child[0]-> nu_parent = 2;
  root->child[1]->color = GREEN;
  root->child[1]-> nu = 1;
  root->child[1]-> nu_parent = 3;
  root->child[0]->child[0]->color = GREEN;
  root->child[0]->child[0]-> nu = 3;
  root->child[0]->child[0]-> nu_parent = 1;
  root->child[0]->child[1]->color = YELLOW;
  root->child[0]->child[1]-> nu = 3;
  root->child[0]->child[1]-> nu_parent = 1;
  root->child[1]->child[0]->color = RED;
  root->child[1]->child[0]-> nu = 3;
  root->child[1]->child[0]-> nu_parent = 1;
  root->child[1]->child[1]->color = RED;
  root->child[1]->child[1]-> nu = 3;
  root->child[1]->child[1]-> nu_parent = 1;

  int count=0;
  assign_processes(root, 20);
  visualize_tree(root,&opts,count++);
  prune_diverged_nodes (root, &opts);
  construct_viable_paths (root);
  print_tree (root);
  visualize_tree(root,&opts,count++);
  choose_viable_paths (root);
  visualize_tree(root,&opts,count++);
  advance_root_node (&root, &opts);
  visualize_tree(root,&opts,count++);
  printf("Calling delete_tree...\n");
  free (root->z_parent); /* Only one needs deallocation */
  delete_tree(root);
  return 0;
}

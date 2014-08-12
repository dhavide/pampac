#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main () {
  /* Constructs long string of nodes to test advance_root_node. */
  options_struct opts = set_options (/*depth=*/5, /*width=*/1);
  PTnode* root = make_tree (&opts);
  /* Assign values to node fields. */
  root->child[0]->color = GREEN;
  root->child[0]->child[0]->color = GREEN;
  root->child[0]->child[0]->child[0]->color = GREEN;
  root->child[0]->child[0]->child[0]->child[0]->color = YELLOW;

  assign_processes(root, 1);
  int count = 0;
  visualize_tree(root,&opts,count++);
  advance_root_node (&root, &opts);
  visualize_tree(root,&opts,count++);
  free (root->z_parent);
  delete_tree(root);
  return 0;
}

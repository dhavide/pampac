#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);
extern void tree43 (PTnode* root);

int
main () {
  options_struct opts = set_options (4,3);
  PTnode* root = make_tree (&opts);
  tree43 (root);
  assign_processes(root, 20);
  int count=0;
  visualize_tree(root,&opts,count++);
  print_tree (root);
  prune_diverged_nodes(root, &opts);
  visualize_tree(root,&opts,count++);
  print_tree (root);
  printf("Calling delete_tree...\n");
  delete_tree(root);
  return 0;
}

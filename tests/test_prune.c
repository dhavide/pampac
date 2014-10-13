#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);
extern void tree43 (PTnode* root);

int
main () {
  options_struct opts = set_options (4,3);
  opts.verbose = 5;
  printf ("%s: opts.input_filename = %s\n", __func__, opts.input_filename);
  PTnode* root = make_tree (&opts);
  tree43 (root);
  assign_processes (root, 20);

  visualize_tree(root,&opts);
//  print_tree (root);
  prune_diverged_nodes(root, &opts);
  visualize_tree(root,&opts);
//  print_tree (root);
//  print_tree (root);
  printf ("Cleaning up...\n");
  delete_tree (root);
  free (opts.scale_factors);
  return 0;
}

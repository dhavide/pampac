#include "pampac.h"
/**********************************************************************/
/* This function advances the root node in the event that it has only */
/* one child node and that node has converged (i.e., is CONVERGED). This  */
/* updating iterates to successive depths in the event that there are */
/* several successive depths with solitary CONVERGED nodes.               */
/**********************************************************************/
void
advance_root_node (PTnode ** root, options_struct * opts) {
  int k;
  bool has_solitary_child, has_converged, has_succeeded;
  PTnode *root_tmp = *root;
  /* Replace root node iff it has a solitary CONVERGED child. */
  has_solitary_child = count_children (root_tmp) == 1;
  has_converged = root_tmp->state == CONVERGED;
  while ( has_solitary_child && has_converged ) {
    /* Identify the index "k" of the single non-NULL child. */
    for (k = 0; k < root_tmp->max_children; k++)
      if (root_tmp->child[k]!=NULL)
        break;
    if (root_tmp->child[k]->state != CONVERGED)
      break;
    /* Dump old root data to disk & reassign the root node. */
    has_succeeded = write_root_coordinates (root_tmp, opts);
    /* Update root_tmp's pointers to become new root node. */
    root_tmp = root_tmp->child[k];
    /* Carefully free pointers associated with old root node. */
    (*root)->child[k] = NULL;
    free ((*root)->z_init);
    (*root)->z_init = NULL;
    free ((*root)->T_init);
    (*root)->T_init = NULL;
    free ((*root)->child);
    (*root)->child = NULL;
    *root = root_tmp;
    if (opts->verbose>0)
      printf ("Advancing root node: lambda=%10g\n",
              (*root)->z[opts->lambda_index]);
    /* Check whether new root node has a solitary CONVERGED child. */
    has_solitary_child = count_children (root_tmp) == 1;
    has_converged = root_tmp->state == CONVERGED;
  }
  assign_depth (*root, 0); /* Update to reflect changes in root node. */
  return;
}

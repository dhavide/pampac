#include "pampac.h"
/**********************************************************************/
/* This function advances the root node in the event that it has only */
/* one child node and that node has converged (i.e., is GREEN). This  */
/* updating iterates to successive depths in the event that there are */
/* several successive depths with solitary GREEN nodes.               */
/**********************************************************************/
void
advance_root_node (PTnode ** root, options_struct * opts)
{
  int k;
  PTnode *tmp_node = *root;
  /* Replace root node if it has a solitary GREEN child. */
  while ((count_children (tmp_node) == 1) && (tmp_node->color==GREEN))
    {
      /* Identify the index "k" of the single non-NULL child. */
      for (k = 0; k < tmp_node->max_children; k++)
	if (tmp_node->child[k]!=NULL)
	  break;
      if (tmp_node->child[k]->color != GREEN)
	break;
      /* Dump old root data to disk & reassign the root node. */
      write_root_coordinates (tmp_node, opts);
      tmp_node = tmp_node->child[k];
      (*root)->child[k] = NULL;
      free ((*root)->z_parent); (*root)->z_parent = NULL;
      free ((*root)->T_parent); (*root)->T_parent = NULL;
      free ((*root)->child);    (*root)->child = NULL;
      *root = tmp_node;
      if (opts->verbose>0)
	printf ("Advancing root node: lambda=%10g\n",
		(*root)->z[opts->lambda_index]);
    }
  assign_depth (*root, 0); /* Update to reflect changes in root node. */
  return;
}

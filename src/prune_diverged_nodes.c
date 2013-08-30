#include "pampac.h"
/**********************************************************************/
/* This function recursively traverses the tree, deleting black nodes */
/* and any associated subtrees.                                       */
/**********************************************************************/
void
prune_diverged_nodes (PTnode *alpha, options_struct *opts)
{
  bool have_rejected_all_children = true;
  bool all_children_NULL = true;
      
  for (int k = 0; k < alpha->max_children; k++)
    {
      PTnode* beta = alpha->child[k];
      if (beta != NULL)
	{
	  all_children_NULL = false;
	  bool is_beta_black = (beta->color == BLACK);
	  have_rejected_all_children &= is_beta_black;
	  if (is_beta_black)
	    {
	      delete_tree (beta);
	      alpha->child[k] = NULL;
	    }
	  else
	    prune_diverged_nodes (beta, opts);
	}
    }
  if (have_rejected_all_children && !all_children_NULL)
    alpha->h *= 0.5;
  return;
}

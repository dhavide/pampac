#include "pampac.h"
/**********************************************************************/
/* This function recursively traverses the tree, deleting black nodes */
/* and any associated subtrees.                                       */
/**********************************************************************/
void
prune_diverged_nodes (PTnode *alpha, options_struct *opts)
{
  int k;
  bool have_pruned_all_children = true;
  bool all_children_NULL = true;
  
  for (k = 0; k < alpha->max_children; k++)
    {
      PTnode* beta = alpha->child[k];
      if (beta != NULL)
	{
	  all_children_NULL = false;
	  bool is_beta_black = (beta->color == BLACK);
	  if (is_beta_black)
	    {
	      delete_tree (beta);
	      alpha->child[k] = NULL;
	    }
	  else
	    prune_diverged_nodes (beta, opts);
	  /* Track whether any children have not been pruned */
	  have_pruned_all_children &= is_beta_black;
	}
    }
  /* If all children were pruned, the parent's step-size h must be
     reduced to ensure all new steps generated are shorter. */
  if (have_pruned_all_children && !all_children_NULL)
    {
      /* Determine minimum and maximum scaling factors */
      double c_min = INFINITY;
      double c_max = -c_min;
      for (k = 0; k < alpha->max_children; k++)
	{
	  double c = opts->scale_process[k];
	  c_min = (c<=c_min)? c : c_min;
	  c_max = (c>=c_max)? c : c_max;
	}
      /* Ensure new step length generates steps that are all 
	 smaller than those that failed. */
      c_min *= (0.9 / c_max);
      if (opts->verbose>1)
	printf ("All children diverged; reduction factor = %g\n", c_min);
      alpha->h *= c_min;
    }
  return;
}

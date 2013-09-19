#include "pampac.h"
/**********************************************************************/
/* Destructor for PTnodes on tree to keep track of computations.      */
/* Recursively traverses tree, freeing nodes and associated memory.   */
/**********************************************************************/
void
delete_tree (PTnode *alpha)
{
  int k;
  if (alpha==NULL)
    return;
  for (k=0; k < alpha->max_children; k++)
    {
      if (alpha->child[k] != NULL)
	{
	  delete_tree (alpha->child[k]);
	}
      alpha->child[k] = NULL;
    }
  free (alpha->child);
  alpha->child = NULL;

  if (alpha->z!=NULL)
    {
      free (alpha->z);
      alpha->z = NULL;
    }
  alpha->z_init = NULL; /* Don't free z_init (shallow copy). */
  if (alpha->T_init!=NULL)
    {
      free (alpha->T_init);
      alpha->T_init = NULL;
    }
  free (alpha);
  alpha = NULL;
}

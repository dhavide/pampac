#include "pampac.h"
/**********************************************************************/
/* Destructor for PTnodes on tree to keep track of computations.      */
/* Recursively traverses tree, freeing nodes and associated memory.   */
/**********************************************************************/
void
delete_tree (PTnode *alpha) {
  int k;
  if (alpha==NULL)
    return;
  if (alpha->child!=NULL) {
    /* First, descend to leaves recursively and delete child nodes */
    for (k=0; k < alpha->max_children; k++) {
      if (alpha->child[k] != NULL)
        delete_tree (alpha->child[k]);
      /* When child nodes are deleted, erase pointers */
      alpha->child[k] = NULL;
    }
    /* Having deleted all child nodes, alpha->child can be released. */
    free (alpha->child);
  }
  /*
   * Now that alpha has no children, release other allocated memory.
   * There are two remaining fields of alpha to release explicitly:
   * the fields z and T_init.
   */
  if (alpha->z!=NULL)
    free (alpha->z);
  if (alpha->T_init!=NULL)
    free (alpha->T_init);
  /*
   * The field z_init points to the field z in alpha's parent node.
   * As such, the memory associated with z_init is released when the
   * parent of alpha is deleted. A special case is the root node; when
   * it is deleted, it is necessary to release the memory associated
   * with the field z_init explicitly at that time.
   */
  if ((alpha->z_init!=NULL) && (alpha->depth==0))
     free (alpha->z_init);
  /*
   * Finally, release memory associated with node alpha itself. Be
   * aware that this does not affect the pointer in the calling stack
   * frame, i.e., "call-by-value" means that the pointer itself will
   * still point to some memory location in that stack frame.
   */
  if (alpha->depth==0)
    printf("delete_tree: Removing root node after deleting tree.\n");
  free (alpha);

  return;
}

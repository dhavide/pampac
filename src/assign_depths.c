#include "pampac.h"
void assign_depths (PTnode *alpha, int depth)
{
  for (int k=0;  k<alpha->max_children; k++)
    {
      PTnode *beta = alpha->child[k];
      if (beta != NULL)
	assign_depths (beta, depth+1 );
    }
  alpha->depth = depth;
}

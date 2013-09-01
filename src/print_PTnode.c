#include "pampac.h"
/**********************************************************************/
/* Given a PTnode, print its contents to stdout.                      */
/**********************************************************************/
void
print_PTnode (PTnode *alpha)
{
  printf ("\nN_dim = %d\n", alpha->N_dim);
  printf ("depth = %d\n", alpha->depth);
  printf ("max_children = %d\n", alpha->max_children);
  printf ("label = %d\n", alpha->label);
  printf ("pid = %d\n", alpha->pid);
  printf ("color = ");
  print_color (alpha,stdout);
  printf("\n");
  printf ("nu = %d\n", alpha->nu);
  printf ("nu_parent = %d\n", alpha->nu_parent);
  printf ("nu_valid = %d\n", alpha->nu_valid);
  printf ("nu_viable = %d\n", alpha->nu_viable);
  printf ("h_init = %g\n", alpha->h_init);
  printf ("h = %g\n", alpha->h);
  printf ("res_norm = %g\n", alpha->res_norm);
  printf ("valid_path_length = %g\n", alpha->valid_path_length);
  printf ("valid_index = %d\n", alpha->valid_index);
  printf ("viable_path_length = %g\n", alpha->viable_path_length);
  printf ("viable_index = %d\n", alpha->viable_index);

  // For sufficiently small problems, print vectors in full
  int k, nt=4;
  if (alpha->z==NULL)
    printf("alpha->z=%p\n",alpha->z);
  else
    if (alpha->N_dim>=nt)
      {
	printf("z = (");
	for (k=0; k<nt; k++)
	  printf(" %10.3e,", alpha->z[k]);
	printf("... )\n");
      }
       
  if (alpha->z_parent==NULL)
    printf("alpha->z_parent=%p\n",alpha->z);
  else
    if (alpha->N_dim>=nt)
      {
	printf("z_parent = (");
	for (k=0; k<nt; k++)
	  printf(" %10.3e,", alpha->z_parent[k]);
	printf("... )\n");
      }

  if (alpha->T_parent==NULL)
    printf("alpha->T_parent=%p\n",alpha->z);
  else
    if (alpha->N_dim>=nt)
      {
	printf("T_parent = (");
	for (k=0; k<nt; k++)
	  printf(" %10.3e,", alpha->T_parent[k]);
	printf("... )\n");
      }
}

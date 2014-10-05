#include "pampac.h"
/**********************************************************************/
/* This routine functions primarily as a wrapper to the user's ouuput */
/* routine write_coordinates. The primary purpose is to pass the data */
/* from the PTnode data structure to the user transparently.          */
/**********************************************************************/
bool
write_root_coordinates (PTnode *node, options_struct *opts) {
  if (opts->verbose>0)
    printf ("Writing point on curve to file: lambda = %12.5g\n",
            node->z[opts->lambda_index]);
  bool has_succeeded = write_coordinates (node->N_dim, node->z);
  return has_succeeded;
}

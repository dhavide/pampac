#include "pampac.h"
/**********************************************************************/
/* Dumps data for a point on the curve into a file.                   */
/**********************************************************************/
void
write_root_coordinates (PTnode *node, options_struct *opts) {
  int k;
  FILE *out_file;
  out_file = fopen (opts->output_filename, "a+");
  if (out_file == NULL) {
    printf ("Error Opening File.\n");
    exit (1);
  }
  if (opts->verbose>0)
    printf ("Writing point on curve to file: lambda = %12.5g\n",
            node->z[opts->lambda_index]);
  for (k = 0; k < node->N_dim; k++)
    fprintf (out_file, " %4.14e", node->z[k]);
  fprintf (out_file, "\n");
  fclose (out_file);
}

#include "pampac.h"
/**********************************************************************/
/* Loads data for initial point on curve from a file.                 */
/**********************************************************************/
bool
load_initial_coordinates (PTnode* node, options_struct* opts)
{
  int N = opts->N_dim, k = 0;
  char *file_name = opts->input_file_name;
  double *z = node->z;
  FILE *input_file;
  input_file = fopen (file_name, "r");
  if (input_file == NULL)
    {
      printf ("load_initial_coordinates: Error Opening File.\n");
      return (false);
    }
  /* Actually parse the input file. */
  for (k = 0; k < N; k++)
    {
      fscanf (input_file, "%lf", z);
      z++;
    }
  fclose (input_file);
  return (true);
}

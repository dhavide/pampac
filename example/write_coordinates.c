#include <stdio.h>
#include <stdbool.h>
/**********************************************************************/
/* This is a user-specified routine to output computed results as the */
/* user chooses. The routine receives an integer (N_dim) and an array */
/* of doubles z (equivalently, a pointer to a double) as input        */
/* arguments. The value returned is a boolean value to indicate       */
/* whether the coordinates were successfully output as desired (this  */
/* flag is used in the principal iteration in the PAMPAC library to   */
/* gracefully halt execution in the event of disk errors, etc. This   */
/* file is intended to serve as a template from which users can write */
/* and tailor outpu as they desire.                                   */
/**********************************************************************/

bool write_coordinates (int N_dim, double *z) {
  bool has_succeeded = true; /* Writing succeeds by default */
  int fprintf_flag;
  FILE *out_file;
  char output_filename[] = "data/output.txt";

  /* We choose in this instance to write the elements of z row-wise
   * into an ASCII file. Successive points computed on the curve are
   * appended as subsequent lines in the file out_file, so the file
   * pointer is opened with the flag "a+". Users may choose to do
   * otherwise as suits their computational needs. */

  out_file = fopen (output_filename, "a+");
  if (out_file == NULL) {
    printf ("write_coordinates: Eror opening file '%s'.\n",
             output_filename);
    return false;
  }

  /* Loop that actually writes to file. Observe the use of the integer 
   * to catch errors in file-writing so that the exception can be
   * caught from calling routines in the PAMPAC library. */
  for (int k = 0; k < N_dim; k++) {
    fprintf_flag = fprintf (out_file, " %4.14e", z[k]);
    if (fprintf_flag<0) {
      printf ("write_coordinates: Aborting, error during file-write\n");
      has_succeeded = false;
      break;
    }
  }
  /* terminate file with a newline (again verifying success). */
  fprintf_flag = fprintf (out_file, "\n");
  if (fprintf_flag<0) {
    printf ("write_coordinates: Aborting, error during file-write\n");
    has_succeeded = false;
  }

  /* Clean up and inform calling routine of success or failure. */
  fclose (out_file);
  return has_succeeded;
}

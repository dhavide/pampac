#include <stdio.h>
#include <string.h>
#include "pampac.h"
/**********************************************************************/
/* Uses input text file to determine problem and parameters required  */
/* by the algorithm to control and tune performance.                  */
/**********************************************************************/
bool
parse_options (int argc, char *argv[], options_struct *opts) {
  FILE *param_file;
  bool has_succeeded = false;

  if (argc<2) {
    printf("parse_options: Parameter_file required command-line input.\n");
    return false;
  }
  param_file = fopen (argv[1], "r");
  if (param_file == NULL) {
    printf ("parse_options: Error Opening File.\n");
    return has_succeeded;
  }

  /* Parse param_file, one line at a time */
  int n_bytes = 0;
  char *line = NULL;
  int bytes_read = getline (&line, &n_bytes, param_file);
  while (bytes_read>0) {
    char param[7];
    bool has_parsed;
    char* parameter_name = malloc (bytes_read*sizeof(char));
    char* value = malloc (bytes_read*sizeof(char));
    if ( sscanf (line, "%s %s %s", param, parameter_name, value)==3 &&
         strcmp (param, "@param") == 0 )
      has_parsed = assign_options (parameter_name, value, opts);
    free (parameter_name);
    parameter_name = NULL;
    free (value);
    value = NULL;
    if (!has_parsed) {
      printf("parse_options: Invalid option parsed.\n");
      fclose (param_file);
      return has_succeeded;
    }
    bytes_read = getline (&line, &n_bytes, param_file);
  }
  fclose (param_file);

  /* Set opts.lambda_dir using opts.h_init */
  double h = opts->h_init;
  /* yields +1 or -1 (direction of initial step along curve) */
  opts->lambda_dir = (h > 0) ? +1.0 : -1.0;
  opts->h_init = (h > 0) ? h : -h;

  has_succeeded = validate_options (opts);
  return has_succeeded;
}



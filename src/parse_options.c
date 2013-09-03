#include <string.h>
#include "pampac.h"
/**********************************************************************/
/* Uses input text file to determine problem and parameters required  */
/* by the algorithm to control and tune performance.                  */
/**********************************************************************/
options_struct
parse_options (char *file_name)
{
  FILE *param_file;
  options_struct opts;
  char buff[BUFSIZ], def[STRING_LEN];
  char param_name[STRING_LEN], param_value[STRING_LEN];
  bool is_max_children_set = false;

  opts.verbose = 0; /* Do not print output by default */
  param_file = fopen (file_name, "r");
  if (param_file == NULL)
    {
      printf ("Error Opening File.\n");
      exit (1);
    }
  /* Parse param_file, one line at a time */
  while (fgets (buff, sizeof (buff), param_file) != NULL)
    {
      if (sscanf (buff, "%s %s %s", def, param_name, param_value) == 3
	  && strcmp (def, "@param") == 0)
	{
	  /* Expect positive integer */
	  if (strcmp (param_name, "N_DIM") == 0)
	    opts.N_dim = atoi (param_value);
	  /* Expect any real value */
	  if (strcmp (param_name, "LAMBDA_MIN") == 0)
	    opts.lambda_min = atof (param_value);
	  /* Expect any real value */
	  if (strcmp (param_name, "LAMBDA_MAX") == 0)
	    opts.lambda_max = atof (param_value);
	  /* Expect integer in range [0,N_dim-1] */
	  if (strcmp (param_name, "LAMBDA_INDEX") == 0)
	    opts.lambda_index = atoi (param_value);
	  /* Expect (small) positive real value */
	  if (strcmp (param_name, "DELTA_LAMBDA") == 0)
	    opts.delta_lambda = atof (param_value);
	  /* Expect (small) positive real value */
	  if (strcmp (param_name, "H_MIN") == 0)
	    opts.h_min = atof (param_value);
	  /* Expect (small) positive real value */
	  if (strcmp (param_name, "H_MAX") == 0)
	    opts.h_max = atof (param_value);
	  /* Expect (small) positive real value */
	  if (strcmp (param_name, "H_INIT") == 0)
	    opts.h_init = atof (param_value);
	  /* Expect positive integer */
	  if (strcmp (param_name, "MAX_ITER") == 0)
	    opts.max_iter = atoi (param_value);
	  /* Expect (small) positive real value */
	  if (strcmp (param_name, "TOL_RESIDUAL") == 0)
	    opts.tol_residual = atof (param_value);
	  /* Expect (moderate) positive real value */
	  if (strcmp (param_name, "GAMMA") == 0)
	    opts.gamma = atof (param_value);
	  /* Expect (moderate) positive real value */
	  if (strcmp (param_name, "MU") == 0)
	    opts.mu = atof (param_value);
	  /* Expect positive integer */
	  if (strcmp (param_name, "MAX_DEPTH") == 0)
	    opts.max_depth = atoi (param_value);
	  /* Expect positive integer */
	  if (strcmp (param_name, "VERBOSE") == 0)
	    opts.verbose = atoi (param_value);
	  /* Expect positive integer */
	  if (strcmp (param_name, "MAX_GLOBAL_ITER") == 0)
	    opts.max_global_iter = atoi (param_value);

	  /* Expect positive integer */
	  if (strcmp (param_name, "MAX_CHILDREN") == 0)
	    {
	      opts.max_children = atoi (param_value);
	      is_max_children_set = true;
	    }

	  if (is_max_children_set)
	    {
	      int i;
	      for (i = 0; i <= opts.max_children; i++)
		{
		  char str[19] = "SCALE_PROCESS_";
		  char str_num[3];
		  sprintf (str_num, "%d", i);
		  strcat (str, str_num);
		  /* Expect positive real value */
		  if (strcmp (param_name, str) == 0)
		    opts.scale_process[i] = atof (param_value);
		}
	    }
	}
    }

  /* Set opts.lambda_dir using opts.h_init */
  /* yields +1 or -1 (direction of initial step along curve) */
  double h = opts.h_init;
  opts.lambda_dir = (h > 0) ? +1.0 : -1.0;
  opts.h_init = abs(h);

  fclose (param_file);
  return opts;
}

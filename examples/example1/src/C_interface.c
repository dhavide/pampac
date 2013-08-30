#include <stdio.h>
#include "C_interface.h"
void single_corrector_step(int N, double *y, double *Tangent)
{
	double out_vars[3];
	int errFlag = FALSE;
  /* Extra underscore in identifier "single_corrector_step_f_" required
     when linking to FORTRAN object code (GNU compilers) */

	single_corrector_step_f_( y, Tangent, out_vars, &errFlag);
   //! Actally, errors should be caught here, too!
}

/***********************************************************************/
void compute_residual(int N, double *y, double *res)
{
  int errflag = FALSE;
  /* Extra underscore in identifier "compute_residual_f_" required
     when linking to FORTRAN object code (GNU compilers) */
  calculate_residual_f_( y, res, &errflag );

  /* Catch error flag: ensure process killed by returning large residual */
  if (errflag == 0)
    {
      printf ("Error in calculate_residual_f_ !");
      return;
    }
  else
    return;
}

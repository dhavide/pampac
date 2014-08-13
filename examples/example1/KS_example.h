#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#include <clapack.h>
#include "pampac.h"

extern void setup_globals (int);
extern void teardown_globals ();
extern void set_differentiation_matrices (int);
extern void compute_Jacobian (int, fftw_complex, fftw_complex, double*);

double Aval, *Dmatrix2, *Dmatrix4, *Res;
fftw_complex *Dmatrix, *X_cplx, *DX_cplx, *RHS_cplx, *Jac_cplx;
int *ipiv;
fftw_plan plan_fft, plan_ifft;

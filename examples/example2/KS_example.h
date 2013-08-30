#include <math.h>
#include <complex.h>
#include <fftw3.h>

extern void print_cplx (int N_cplx, fftw_complex * z);
extern void print_real (int N_real, double *z);
extern void set_differentiation_matrices (int N_grid, fftw_complex *D, double *D2, double *D4);
extern void setup_teardown (int N_real);

double Aval;

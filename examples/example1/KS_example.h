#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "pampac.h"

extern void print_cplx (int, fftw_complex*);
extern void print_real (int, double*);
extern void set_differentiation_matrices (int, fftw_complex*, double*, double*);
extern void setup_teardown (int);

double Aval;

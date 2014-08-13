#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#include <clapack.h>
#include "pampac.h"

extern void set_differentiation_matrices (int, fftw_complex*, double*, double*);
extern void setup_teardown (int);

double Aval;

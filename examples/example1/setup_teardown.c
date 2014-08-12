#include "KS_example.h"

void setup_teardown(int N) {
  // Required to create or destroy FFT plans for both routines
  compute_residual(N, NULL, NULL);
  single_corrector_step(N, NULL, NULL);
  return;
}

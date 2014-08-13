#include "KS_example.h"
/**********************************************************************/
void setup_globals(int N_dim) {
  /* N_dim = number of REAL variables input */
  /* N_grid = number of collocation points  */
  int N_grid = (N_dim / 2) - 2;
  /* Setting up static arrays and data structures */
  Dmatrix = fftw_alloc_complex (N_grid);
  Dmatrix2 = fftw_alloc_real (N_grid);
  Dmatrix4 = fftw_alloc_real (N_grid);
  set_differentiation_matrices (N_grid);
  Res = fftw_alloc_real (N_dim-1);
  X_cplx = fftw_alloc_complex (N_grid);
  DX_cplx = fftw_alloc_complex (N_grid);
  Jac_cplx = fftw_alloc_complex (pow(N_grid+2,2));
  RHS_cplx = fftw_alloc_complex (N_grid+2);
  ipiv = malloc ((N_grid+2) * sizeof(*ipiv));
  plan_fft = fftw_plan_dft_1d (N_grid, X_cplx, X_cplx,
                               FFTW_FORWARD, FFTW_MEASURE);
  plan_ifft = fftw_plan_dft_1d (N_grid, X_cplx, X_cplx,
                                FFTW_BACKWARD, FFTW_MEASURE);
  return;
}
/**********************************************************************/
void teardown_globals() {
  /* Cleaning up static arrays and data structures */
  fftw_free (Dmatrix);
  Dmatrix = NULL;
  fftw_free (Dmatrix2);
  Dmatrix2 = NULL;
  fftw_free (Dmatrix4);
  Dmatrix4 = NULL;
  fftw_free (Res);
  Res = NULL;
  fftw_free (X_cplx);
  X_cplx = NULL;
  fftw_free (DX_cplx);
  DX_cplx = NULL;
  fftw_free (Jac_cplx);
  Jac_cplx = NULL;
  fftw_free (RHS_cplx);
  RHS_cplx = NULL;
  free (ipiv);
  ipiv = NULL;
  fftw_destroy_plan (plan_fft);
  fftw_destroy_plan (plan_ifft);
  return;
}
/**********************************************************************/
void
set_differentiation_matrices (int N_grid)
/**********************************************************************/
/* Helper routine for compute_residual and for single_corrector_step. */
/* Observe that differentiation matrices are diagonal and hence are   */
/* stored as 1D arrays. Also notice that the array Dmatrix is purely  */
/* imaginary, so Dmatrix2==(Dmatrix)^2 and Dmatrix4==(Dmatrix)^4 are  */
/* both real matrices (and hence are declared as doubles).            */
/* Note: it is assumed that memory allocation occurs elsewhere.       */
/**********************************************************************/
{
  int k;
  Dmatrix[0] = 0.0;
  Dmatrix2[0] = 0.0;
  Dmatrix4[0] = 0.0;
  for (k = 1; k < N_grid / 2; k++) {
    Dmatrix[k] = 1.0I*k;
    Dmatrix[N_grid - k] = -Dmatrix[k];
    Dmatrix2[k] = -pow (k, 2);
    Dmatrix2[N_grid - k] = -pow (k, 2);
    Dmatrix4[k] = pow (k, 4);
    Dmatrix4[N_grid - k] = pow (k, 4);
  }

  Dmatrix[N_grid / 2] = 0.0;
  Dmatrix2[N_grid / 2] = -pow (N_grid, 2) / 4.0;
  Dmatrix4[N_grid / 2] = pow (N_grid, 4) / 16.0;
  return;
}
/**********************************************************************/
void compute_Jacobian (int N_dim, fftw_complex c,
                         fftw_complex nu, double * T) {
  int N_grid, k, ell, index;
  /* N_dim = number of REAL variables input */
  /* N_grid = number of collocation points  */
  N_grid = (N_dim / 2) - 2;
  /* Resetting Jacobian entries to 0.0 */
  for (k=0; k<N_grid+2; k++)
    for (ell=0; ell<N_grid+2; ell++)
      Jac_cplx[k+(N_grid+2)*ell] = 0.0;

  for (k=0; k<N_grid; k++) {
    /* Setting diagonal elements */
    index = k + (N_grid+2) * k;
    Jac_cplx[index] = (-c) * Dmatrix[k];
    Jac_cplx[index] += Dmatrix2[k];
    Jac_cplx[index] += nu*Dmatrix4[k];
    /* Derivative wrt wavespeed c */
    index = k + (N_grid+2) * N_grid;
    Jac_cplx[index] = -Dmatrix[k] * X_cplx[k];
    /* Derivative wrt viscosity nu */
    index = k + (N_grid+2) * (N_grid+1);
    Jac_cplx[index] = Dmatrix4[k] * X_cplx[k];
    /* 2nd to last row = (D*X)' */
    index = N_grid + (N_grid+2) * k;
    Jac_cplx[index] = conj( Dmatrix[k] * X_cplx[k] );
    /* last row = T' */
    index = (N_grid+1) + (N_grid+2) * k;
    Jac_cplx[index] = T[2*k] - 1.I*T[2*k+1];
  }

  /* Complete last two entries of last row */
  index = (N_grid+1) + (N_grid+2) * N_grid;
  Jac_cplx[index] = T[N_dim-4]-1.0I*T[N_dim-3];
  index = (N_grid+1) + (N_grid+2) * (N_grid+1);
  Jac_cplx[index] = T[N_dim-2]-1.0I*T[N_dim-1];

  /* Accumulate advective terms in Jacobian matrix */
  for (k=0; k<=N_grid-1; k++)
    for (ell=0; ell<=k; ell++) {
        index = k + (N_grid+2) * ell;
      Jac_cplx[index] += Dmatrix[k-ell] * X_cplx[k-ell] / N_grid;
      Jac_cplx[index] += Dmatrix[ell] * X_cplx[k-ell] / N_grid;
    }

  for (k=0; k<=N_grid-2; k++)
    for (ell=k+1; ell<=N_grid-1; ell++) {
        index = k + (N_grid+2) * ell;
      Jac_cplx[index] += Dmatrix[k-ell+N_grid] *
                          X_cplx[k-ell+N_grid] / N_grid;
      Jac_cplx[index] += Dmatrix[ell] * X_cplx[k-ell+N_grid] / N_grid;
    }

  /* Use DFTs to compute nonlinear terms of function. */
  /* Equivalent to X_cplx = N_grid * ifft(X_cplx) in Matlab... */
  fftw_execute_dft (plan_ifft, X_cplx, X_cplx);
  for (k=0; k<N_grid; k++) {
    /* Rescale: FFTW's "ifft" is scaled by N_grid */
    X_cplx[k] /= N_grid;
    /* Same as "X_cplx = Aval*cos(ifft(X_cplx))" in Matlab */
    X_cplx[k] = Aval*ccos( X_cplx[k] );
  }
  /* Equivalent to X_cplx = fft(X_cplx) in Matlab... */
  fftw_execute_dft (plan_fft, X_cplx, X_cplx);

  /* Accumulate nonlinear trigonometric derivative terms in Jacobian matrix */
  for (k=0; k<=N_grid-1; k++)
    for (ell=0; ell<=k; ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell] / N_grid;

  for (k=0; k<=N_grid-2; k++)
    for (ell=k+1; ell<=N_grid-1; ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell+N_grid] / N_grid;
return;
}

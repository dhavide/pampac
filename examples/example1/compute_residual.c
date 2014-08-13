#include "KS_example.h"
/******************************************************************************/

void
set_differentiation_matrices (int N_grid,
                              fftw_complex * D, double *D2, double *D4)
/******************************************************************************/
/* Helper routine for compute_residual and for single_corrector_step. Observe */
/* that differentiation matrices are diagonal and so are stored as 1D arrays. */
/* Also notice that the array D is pure imaginary, so D2==D^2 and D4==D^4 are */
/* both real matrices (and hence are declared as doubles).                    */
/* Note: it is assumed that memory allocation ocurs in the calling workspace. */
/******************************************************************************/
{
  int k;
  D[0] = 0.0;
  D2[0] = 0.0;
  D4[0] = 0.0;
  for (k = 1; k < N_grid / 2; k++) {
    D[k] = 1.0I*k;
    D[N_grid - k] = -D[k];
    D2[k] = -pow (k, 2);
    D2[N_grid - k] = -pow (k, 2);
    D4[k] = pow (k, 4);
    D4[N_grid - k] = pow (k, 4);
  }

  D[N_grid / 2] = 0.0;
  D2[N_grid / 2] = -pow (N_grid, 2) / 4.0;
  D4[N_grid / 2] = pow (N_grid, 4) / 16.0;
  return;
}

/******************************************************************************/

void
compute_residual (int N_real, const double *Z, double *Res)
/******************************************************************************/
/* This function computes a particular nonlinear function of Z and returns    */
/* the result in the array Res. The function arises from a pseudo-spectral    */
/* discretisation of the one-dimensional ODE                                  */
/*                                                                            */
/* -c \psi' + \psi \psi' + \psi'' + \nu\psi^{(iv)} + A\sin(\psi) = 0,         */
/* \psi(0) = \psi(2\pi).                                                      */
/*                                                                            */
/* The BVP is discretized by looking for trigonometric polynomial solutions   */
/* \psi(\xi) = \sum_{k=0}^{N_{grid}-1} X_{k} \exp(i k \xi)                    */
/* and using collocation on a uniform grid to determine the coefficients X_k. */
/* That is, the collocation points are taken to be \xi_\ell=2\pi\ell/N_{grid}.*/
/* The resulting nonlinear vector field can be evaluated efficiently using    */
/* discrete Fourier transforms (hence the use of the FFTW library).           */
/*                                                                            */
/* This BVP arises in seeking travelling wave solutions of the form psi(x-ct) */
/* for the time-dependent PDE                                                 */
/*                                                                            */
/* u_t + u u_x + u_{xx} + \nu u_{xxxx} + A \sin(u) = 0                        */
/*                                                                            */
/* (a modified Kuramoto-Shivashinsky equation) with periodic BCs              */
/*                                                                            */
/* u(0,t) = u(2*\pi,t) for all t.                                             */
/*                                                                            */
/* The input vector Z consists of N_{real}=2 N_{grid} + 4 components: in      */
/* order, the components of Z are the real and imaginary parts of X_{k}, the  */
/* Fourier coefficients in alternating order, the wave speed c (real and      */
/* imaginary parts) and the viscosity \nu (real and imaginary parts).         */
/* To clarify, the components of Z are as follows:                            */
/*                                                                            */
/* Z[0] = Re(X[0]) <- 0th Fourier coefficient                                 */
/* Z[1] = Im(X[0]) <- 0th Fourier coefficient                                 */
/* Z[2] = Re(X[1]) <- 1st Fourier coefficient                                 */
/* Z[3] = Im(X[1]) <- 1st Fourier coefficient                                 */
/*   ...                                                                      */
/* Z[2*N_grid-2] = Re(X[N_grid-1]) <- (N_grid-1)st Fourier coefficient        */
/* Z[2*N_grid-1] = Im(X[N_grid-1]) <- (N_grid-1)st Fourier coefficient        */
/* Z[2*N_grid]   = Re(c)   <- wavespeed                                       */
/* Z[2*N_grid+1] = Im(c)   <- wavespeed                                       */
/* Z[2*N_grid+2] = Re(nu)  <- viscosity                                       */
/* Z[2*N_grid+3] = Im(nu)  <- viscosity                                       */
/*                                                                            */
/* The computed vector field is returned in the array of doubles Res that has */
/* 2*N_{grid}+3 components. The array Res also uses a compressed complex      */
/* storage scheme (i.e., alternating real and imaginary parts) to represent   */
/* a complex vector using an array of doubles. The last three components of   */
/* are set to zero.                                                           */
/*                                                                            */
/* Passing NULL pointers for either Res or Z signals to set up or clean up    */
/* static data structures used internally.                                    */
/******************************************************************************/
{
  static fftw_complex *D = NULL;
  static double *D2 = NULL, *D4 = NULL;
  static fftw_complex *X_cplx = NULL, *DX_cplx = NULL, *Res_cplx = NULL;
  static fftw_plan plan_fft, plan_ifft;
  static bool setup_complete = false;
  int N_grid, k;
  fftw_complex c, nu;

  /* N_real = number of REAL variables input */
  N_grid = (N_real / 2) - 2;    /* Number of collocation points */
  /*
     A NULL pointer for the input array Res initiates prelimary set up or
     final clean up. Static memory is allocated for intermediate work and
     the FFTW plans are initialized. The boolean setup_complete distinguishes
     whether to initialize or release memory & data structures.
   */

  if (Z == NULL || Res == NULL) {
    if (!setup_complete) {
      /* Setting up static arrays and data structures */
      D = fftw_alloc_complex (N_grid);
      D2 = fftw_alloc_real (N_grid);
      D4 = fftw_alloc_real (N_grid);
      X_cplx = fftw_alloc_complex (N_grid);
      DX_cplx = fftw_alloc_complex (N_grid);
      Res_cplx = fftw_alloc_complex (N_grid);

      set_differentiation_matrices (N_grid, D, D2, D4);
      plan_fft =
        fftw_plan_dft_1d (N_grid, X_cplx, X_cplx, FFTW_FORWARD,
                          FFTW_MEASURE);
      plan_ifft =
        fftw_plan_dft_1d (N_grid, X_cplx, X_cplx, FFTW_BACKWARD,
                          FFTW_MEASURE);

      /* Static flag changed to signal future calls to compute_residual */
      setup_complete = true;

      return;
    }
    /* Cleaning up static arrays and data structures */
    fftw_free (D);
    D = NULL;
    fftw_free (D2);
    D2 = NULL;
    fftw_free (D4);
    D4 = NULL;
    fftw_free (X_cplx);
    X_cplx = NULL;
    fftw_free (DX_cplx);
    DX_cplx = NULL;
    fftw_free (Res_cplx);
    Res_cplx = NULL;
    fftw_destroy_plan (plan_fft);
    fftw_destroy_plan (plan_ifft);
    /* Static flag reset for future calls to compute_residual */
    setup_complete = false;
    return;
  }
  /* Main routine: compute nonlinear residual */
  if (!setup_complete) {
    printf
    ("Error: Call made to compute_residual without prior initialisation!\n");
    printf
    ("       Call compute_residual ( N, Z, Res ) with NULL pointer for Z or Res to initialise.\n");
    return;
  }
  c = Z[N_real - 4] + 1.0I*Z[N_real-3];
  nu = Z[N_real - 2] + 1.0I*Z[N_real-1];
  /* Accumulate linear terms first */
  for (k = 0; k < N_grid; k++) {
    X_cplx[k] = Z[2*k] + 1.0I * Z[2*k + 1];
    Res_cplx[k] = D2[k] + nu * D4[k];
    Res_cplx[k] *= X_cplx[k];
    DX_cplx[k] = D[k] * X_cplx[k];
    Res_cplx[k] -= c * DX_cplx[k];
  }
  /* Now use DFTs to compute nonlinear terms of function. */
  /* Equivalent to X_cplx = N_grid * ifft(X_cplx) in Matlab... */
  fftw_execute_dft (plan_ifft, X_cplx, X_cplx);
  /* Equivalent to DX_cplx = N_grid * ifft(DX_cplx) in Matlab... */
  fftw_execute_dft (plan_ifft, DX_cplx, DX_cplx);

  for (k = 0; k < N_grid; k++) {
    X_cplx[k] /= N_grid;    /* Rescale: FFTW's "ifft" is scaled by N_grid */
    DX_cplx[k] /= N_grid;   /* Rescale: FFTW's "ifft" is scaled by N_grid */
    DX_cplx[k] *= X_cplx[k];    /* Same as DX_cplx = ifft(DX_cplx) .* ifft(X_cplx) */
    X_cplx[k] = csin (X_cplx[k]);   /* Same as X_cplx = sin( ifft(X_cplx) ) */
  }

  /* Equivalent to X_cplx = fft(X_cplx) in Matlab */
  fftw_execute_dft (plan_fft, X_cplx, X_cplx);
  /* Equivalent to DX_cplx = fft(DX_cplx) in Matlab */
  fftw_execute_dft (plan_fft, DX_cplx, DX_cplx);

  /* Finally, accumulate nonlinear terms into residual vector and unpack real
     and imaginary components into vector Res of *doubles*. */
  for (k = 0; k < N_grid; k++) {
    Res_cplx[k] += DX_cplx[k];
    Res_cplx[k] += Aval * X_cplx[k];
    Res[2*k] = creal (Res_cplx[k]);
    Res[2*k + 1] = cimag (Res_cplx[k]);
  }
  Res[N_real-4] = 0.0;
  Res[N_real-3] = 0.0;
  Res[N_real-2] = 0.0;
  return;
}
/******************************************************************************/

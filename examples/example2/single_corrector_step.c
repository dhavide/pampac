#include <stdbool.h>
#include <stdlib.h>
#include <clapack.h>
#include "KS_example.h"
void
single_corrector_step (int N_real, double *Z, double *T)
/******************************************************************************/
/* This function implements a single Newton corrector step for the BVP        */
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
/* Z[2*N_grid)   = Re(c)   <- wavespeed                                       */
/* Z[2*N_grid+1] = Im(c)   <- wavespeed                                       */
/* Z[2*N_grid+2] = Re(nu)  <- viscosity                                       */
/* Z[2*N_grid+3] = Im(nu)  <- viscosity                                       */
/*                                                                            */
/* The input array T represents the tangent vector to the curve Z lies on.    */
/* The array of doubles T is also ordered using compressed complex storage.   */
/* The result of the computation modifies the array Z in place, i.e., Z  is   */
/* overwritten upon return from this function.                                */
/*                                                                            */
/* Passing NULL pointers for either Z or T signals to set up or clean up      */
/* static data structures used internally.                                    */
/*                                                                            */
/******************************************************************************/
{
  static fftw_complex *D = NULL;
  static double *D2 = NULL, *D4 = NULL, *Res = NULL;
  static fftw_complex *X_cplx = NULL, *RHS_cplx = NULL, *Jac_cplx = NULL;
  static fftw_plan plan_fft, plan_ifft;
  static int *ipiv=NULL;
  static bool setup_complete = false;
  int N_grid, k, ell, info;
  fftw_complex c, nu;
  /* N_real = number of REAL variables input */

  N_grid = (N_real / 2) - 2;	// Number of collocation points
  /*
     Receiving a NULL pointer for Z or T initiates prelimary set up or
     final clean up. Static memory is allocated for intermediate work and
     the FFTW plans are initialized. The boolean setup_complete distinguishes
     whether to initialize or release memory & data structures.
   */
  if (Z == NULL || T == NULL)
    {
      if (!setup_complete)
	{
	  /* Setting up static arrays and data structures */
	  D = fftw_alloc_complex (N_grid);
	  D2 = fftw_alloc_real (N_grid);
	  D4 = fftw_alloc_real (N_grid);
	  X_cplx = fftw_alloc_complex (N_grid);
	  Res = fftw_alloc_real (2*N_grid);
	  RHS_cplx = fftw_alloc_complex (N_grid+2);
	  Jac_cplx = fftw_alloc_complex (pow(N_grid+2,2));
	  ipiv = malloc((N_grid+2)*sizeof(*ipiv));
	  set_differentiation_matrices (N_grid, D, D2, D4);
	  plan_fft =
	    fftw_plan_dft_1d (N_grid, X_cplx, X_cplx, FFTW_FORWARD,
			      FFTW_MEASURE);
	  plan_ifft =
	    fftw_plan_dft_1d (N_grid, X_cplx, X_cplx, FFTW_BACKWARD,
			      FFTW_MEASURE);
	  /* Static flag changed for future calls to single_corrector_step */
	  setup_complete = true;
	  return;
	}
      /* Cleaning up static arrays and data structures */
      fftw_destroy_plan (plan_fft);
      fftw_destroy_plan (plan_ifft);
      fftw_free (D);
      D = NULL;
      fftw_free (D2);
      D2 = NULL;
      fftw_free (D4);
      D4 = NULL;
      fftw_free (Res);
      Res = NULL;
      fftw_free (X_cplx);
      X_cplx = NULL;
      fftw_free (Jac_cplx);
      Jac_cplx = NULL;
      fftw_free (RHS_cplx);
      RHS_cplx = NULL;
      free (ipiv);
      ipiv = NULL;
      /* Reset static flag for future calls to single_corrector_step */
      setup_complete = false;

      return;
    }
  /* Main routine: compute corrector step */
  if (!setup_complete)
    {
      printf
	("Error: Call made to single_corrector_step without prior initialisation!\n");
      printf
	("       Call single_corrector_step ( N, Z, T ) with NULL pointer for Z or T to initialise.\n");
      return;
    }
  for (k = 0; k < N_grid; k++)
      X_cplx[k] = Z[2*k] + 1.0I * Z[2*k + 1];
  c = Z[N_real - 4] + 1.0I*Z[N_real-3];
  nu = Z[N_real - 2] + 1.0I*Z[N_real-1];
  
  /* Resetting Jacobian entries to 0.0 */
  for (k=0; k<N_grid+2; k++)
    for (ell=0;ell<N_grid+2; ell++)
		Jac_cplx[k+(N_grid+2)*ell] = 0.0;
    
  for (k=0; k<N_grid; k++)
  {
    /* Setting diagonal elements */
    Jac_cplx[k+(N_grid+2)*k] = -c*D[k] + D2[k] + nu*D4[k];
    /* Derivative wrt wavespeed c */
    Jac_cplx[k+(N_grid+2)*N_grid] = -D[k] * X_cplx[k];
    /* Derivative wrt viscosity nu */
    Jac_cplx[k+(N_grid+2)*(N_grid+1)] = D4[k] * X_cplx[k];
    /* 2nd to last row = (D*X)' */
    Jac_cplx[N_grid+(N_grid+2)*k] = conj( D[k]*X_cplx[k] );
    /* last row = T' */
    Jac_cplx[N_grid+1+(N_grid+2)*k] = T[2*k] - 1.I*T[2*k+1];
  }

  /* Complete last two entries of last row */
  Jac_cplx[N_grid+1+(N_grid+2)*N_grid] = T[N_real-4]-1.0I*T[N_real-3];
  Jac_cplx[N_grid+1+(N_grid+2)*(N_grid+1)] = T[N_real-2]-1.0I*T[N_real-1];

  /* Accumulate advective terms in Jacobian matrix */
  for (k=0;k<=N_grid-1;k++)
    for (ell=0;ell<=k;ell++)
      {
	Jac_cplx[k+(N_grid+2)*ell] += D[k-ell] * X_cplx[k-ell] / N_grid;
	Jac_cplx[k+(N_grid+2)*ell] += D[ell] * X_cplx[k-ell] / N_grid;
      }
  
  for (k=0;k<=N_grid-2;k++)
    for (ell=k+1;ell<=N_grid-1;ell++)
      {
	Jac_cplx[k+(N_grid+2)*ell] += D[k-ell+N_grid] * X_cplx[k-ell+N_grid] / N_grid;
	Jac_cplx[k+(N_grid+2)*ell] += D[ell] * X_cplx[k-ell+N_grid] / N_grid;
      }

  /* Use DFTs to compute nonlinear terms of function. */
  /* Equivalent to X_cplx = N_grid * ifft(X_cplx) in Matlab... */
  fftw_execute_dft (plan_ifft, X_cplx, X_cplx);
  for (k=0;k<N_grid;k++)
    {
      /* Rescale: FFTW's "ifft" is scaled by N_grid */
      X_cplx[k] /= N_grid;
      /* Same as "X_cplx = A*cos(ifft(X_cplx))" in Matlab */
      X_cplx[k] = Aval*ccos( X_cplx[k] );
    }
  /* Equivalent to X_cplx = fft(X_cplx) in Matlab... */
  fftw_execute_dft (plan_fft, X_cplx, X_cplx);

  /* Accumulate nonlinear trigonometric derivative terms in Jacobian matrix */
  for (k=0;k<=N_grid-1;k++)
    for (ell=0;ell<=k;ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell] / N_grid;

  for (k=0;k<=N_grid-2;k++)
    for (ell=k+1;ell<=N_grid-1;ell++)
      Jac_cplx[k+(N_grid+2)*ell] += X_cplx[k-ell+N_grid] / N_grid;

  /* Compute the right-hand side and pack into a complex vector */
  compute_residual ( N_real, Z, Res );
  for (k=0;k<N_grid;k++)
    RHS_cplx[k] = Res[2*k] + 1.0I * Res[2*k+1];
  RHS_cplx[N_grid] = 0.0;
  RHS_cplx[N_grid+1] = 0.0;

  /* Solve linear system Jac_cplx*X=RHS_cplx; solution overwrites RHS_cplx */
  info = clapack_zgesv( CblasColMajor, N_grid+2, 1, Jac_cplx, N_grid+2, 
			ipiv, RHS_cplx, N_grid+2);
  if (info  != 0) fprintf(stderr,"failure with error %d\n", info);

  /* Filter out the part of the Newton update that 
     gives complex-valued solutions or parameters
     due to accumulated numerical noise */
  RHS_cplx[0] = creal(RHS_cplx[0]);
  RHS_cplx[N_grid/2] = creal(RHS_cplx[N_grid/2]);
  RHS_cplx[N_grid] = creal(RHS_cplx[N_grid]);
  RHS_cplx[N_grid+1] = creal(RHS_cplx[N_grid+1]);
  for (k=1;k<N_grid/2;k++)
   {
     RHS_cplx[k] = 0.5*(RHS_cplx[k]+conj(RHS_cplx[N_grid-k]));
     RHS_cplx[N_grid-k] = conj(RHS_cplx[k]);
   };


  /* Overwrite array Z using corrector step in; unwrap real & imaginary parts
     (notice minus sign) */
  for (k=0;k<N_grid+2;k++)
    {
      Z[2*k] -= creal(RHS_cplx[k]);
      Z[2*k+1] -= cimag(RHS_cplx[k]);
    }
  return;
}

#include "mpi.h"
#include "string.h"
#include "KS_example.h"
/**********************************************************************/
/* Generic MPI main routine with input options communicated through a */
/* custom data structure opts.                                        */
/**********************************************************************/
int
main (int argc, char *argv[]) {
  int p_id;         /* rank of process */
  int N_p;          /* number of processes */
  int N_grid;       /* number of grid points */
  int N_dim;        /* Number of variables */

  /* Global variable Aval declared in KS_example.h; clumsy, but this
     allows it to be set in one place for use by both functions
     compute_residual and single_corrector_step. */
  N_grid = 2048;
  Aval = -8.089150779042718; /* for N=2048 test case */
  /* Aval = -9.135914993072223;  for N=1024 test case */

  N_dim = 2*(N_grid+2);
  MPI_Init (NULL, NULL);               /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id); /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);  /* find # of processes */
  if (p_id == 0) {
    MPI_Bcast(&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
    setup_globals (N_dim); /* Initialisation */
    master_process (N_p, argc, argv);
    teardown_globals ();   /* Clean-up */

    printf ("PID %d: Cleaning up master process.\n", p_id);
    fflush (stdout);
    MPI_Barrier (MPI_COMM_WORLD);
  } else {
    /* N_dim is read from disk by master and broadcast to slaves */
    MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);

    setup_globals (N_dim); /* Initialisation */
    slave_process (N_dim);
    teardown_globals ();   /* Clean-up */

    printf ("PID %d: Cleaning up slave process.\n", p_id);
    fflush (stdout);
    MPI_Barrier (MPI_COMM_WORLD);
  }
  MPI_Finalize ();      /* shut down MPI */
  return 0;
}

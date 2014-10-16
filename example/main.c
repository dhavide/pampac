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
  int N_dim;        /* number of variables */
  bool success;     /* flag to check for successful setup */

  N_grid = 2048;
  N_dim = 2*(N_grid+2); /* Inelegant but needed */

  MPI_Init (NULL, NULL);                 /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id); /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);  /* find # of processes */
  if (p_id == 0) {
    MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
    success = setup_globals (N_dim); /* Initialisation */
    if (success) {
      master_process (N_p, argc, argv);
      teardown_globals ();   /* Clean-up */
      printf ("PID %d: Cleaning up master process.\n", p_id);
      fflush (stdout);
    }
    MPI_Barrier (MPI_COMM_WORLD); /* Synchronise before shutting down */
  } else {
    /* N_dim is read from disk by master and broadcast to slaves */
    MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
    success = setup_globals (N_dim); /* Initialisation */
    if (success) {
      slave_process (N_dim);
      teardown_globals ();   /* Clean-up */
      printf ("PID %d: Cleaning up slave process.\n", p_id);
      fflush (stdout);
    }
    MPI_Barrier (MPI_COMM_WORLD); /* Synchronise before shutting down */
  }
  MPI_Finalize ();      /* shut down MPI */
  return 0;
}

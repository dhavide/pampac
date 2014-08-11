#include "pampac.h"
#include "mpi.h"
#include "string.h"
#include "KS_example.h"
/**********************************************************************/
/* Generic MPI main routine with input options communicated through a */
/* custom data structure opts.                                        */
/**********************************************************************/
int
main (int argc, char *argv[])
{
  int p_id;         /* rank of process */
  int N_p;          /* number of processes */
  int N_dim;                    /* number of variables */

  /* Global variable Aval declared in KS_example.h; clumsy, but this
     allows it to be set in one place for use by both functions
     compute_residual and single_corrector_step. */
  Aval = -8.089150779042718; /* for N=2048 test case */
  /* Aval = -9.135914993072223;  for N=1024 test case */

  MPI_Init (&argc, &argv);               /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id); /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);  /* find # of processes */
  if (p_id == 0)
    {
      options_struct opts;      /* problem/algorithm parameters */
      opts = parse_options (argv[1]);
      N_dim = opts.N_dim;
      MPI_Bcast(&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
      opts.input_file_name = argv[2];
      opts.output_file_name = argv[3];
      opts.base_name_tree = argv[4];
      /* Initialisation */
      setup_teardown(N_dim);
      master_process (p_id, N_p, &opts);
      setup_teardown(N_dim); /* Clean-up */
      printf("PID %d: Cleaning up master process.\n", p_id);
      fflush(stdout);
      MPI_Barrier (MPI_COMM_WORLD);
    }
  else
    {
      /* N_dim is read from disk by master and broadcast to slaves */
      MPI_Bcast(&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
      setup_teardown(N_dim); /* Initialisation */
      slave_process (N_dim);
      setup_teardown(N_dim); /* Clean-up */
      printf("PID %d: Cleaning up slave process.\n", p_id);
      fflush(stdout);
      MPI_Barrier (MPI_COMM_WORLD);
    }
  MPI_Finalize ();      /* shut down MPI */
  return 0;
}

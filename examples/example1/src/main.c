#include "pampac.h"
#include "mpi.h"
/**********************************************************************/
/* Generic MPI main routine with input options communicated through a */
/* custom data structure opts.                                        */
/**********************************************************************/
int
main (int argc, char *argv[])
{
  int p_id;			/* rank of process */
  int N_p;			/* number of processes */
  int N_dim;                    /* number of variables */

  MPI_Init (&argc, &argv);               /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id); /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);	 /* find # of processes */
  if (p_id == 0)
    {
      printf ("main: p_id = %d\n", p_id);
      options_struct opts;		/* problem/algorithm parameters */
      opts = parse_options (argv[1]);
      N_dim = opts.N_dim;
      MPI_Bcast(&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
      opts.input_file_name = argv[2];
      opts.output_file_name = argv[3];
      opts.base_name_tree = argv[4];
      /* Initialisation */

      master_process (p_id, N_p, &opts);

      printf("PID %d: Cleaning up master process.\n", p_id);
      fflush(stdout);
      MPI_Barrier (MPI_COMM_WORLD);
    }
  else
    {
      /* N_dim is read from disk by master and broadcast to slaves */
      MPI_Bcast(&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);

      slave_process (N_dim);

      printf("PID %d: Cleaning up slave process.\n", p_id);
      fflush(stdout);
      MPI_Barrier (MPI_COMM_WORLD);
    }
  MPI_Finalize ();		/* shut down MPI */
  return 0;
}

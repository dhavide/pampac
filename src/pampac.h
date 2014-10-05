#define EXIT_TAG 202
#define CONTINUE_TAG 201
#ifndef max
#define max( a, b ) (((a)>(b)) ? (a) : (b))
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
/**********************************************************************/
/* Required for extracting required algorithm parameters from a file  */
/**********************************************************************/
typedef struct options_struct {
  int N_dim;                /* dimensions of problem                  */
  double lambda_min;        /* lower continuation parameter value     */
  double lambda_max;        /* upper continuation parameter value     */
  int lambda_index;         /* continuation parameter index [0-based] */
  int lambda_dir;           /* initial direction of tangent (+/- 1)   */
  double delta_lambda;      /* initial perturbation in parameter      */

  double h_min;             /* lower bound on step-size               */
  double h_max;             /* upper bound on step-size               */
  double h_init;            /* starting value for step-size           */

  int max_iter;             /* max. iterations within corrector steps */
  double tol_residual;      /* corrector step convergence tolerance   */
  double mu;                /* threshold for divergence               */
  double gamma;             /* tolerance for "nearing convergence"    */

  int max_depth;            /* limits parallel processes spawned      */
  int max_children;         /* limits chold nodes for processes       */
  double *scale_factors;    /* step-size scaling factors              */
  int max_global_iter;      /* limits global continuation loop        */

  int verbose;              /* integer flag for controlling output    */
  char* input_filename;     /* parameter file from user               */
  char* tree_base_filename; /* used to name tree output files         */
  int tree_filename_num;    /* counter for tree dot files generated   */
} options_struct;

/**********************************************************************/
/* Colors used for distinguishing convergence status of iterations    */
/* running concurrently (represented by nodes of a rooted tree        */
/**********************************************************************/
typedef enum
{ BLACK = 0, GREEN = 1, YELLOW = 2, RED = 3 } NodeColors;
/**********************************************************************/

/**********************************************************************/
/* Computations run concurrently on distinct processors represented   */
/* using PTnode data structure (nodes on a rooted tree)               */
/**********************************************************************/
typedef struct PTnode {
  int N_dim;
  int label;
  int pid;
  int depth;
  NodeColors color;
  int nu;
  int nu_init;
  int nu_valid;
  int nu_viable;
  double h_init;
  double h;
  double res_norm;
  double *z;
  double *T_init;
  double *z_init;
  double valid_path_length;
  double viable_path_length;
  int viable_index;
  int valid_index;
  int max_children;
  struct PTnode **child;
} PTnode;

/**********************************************************************/
/* Queue data structure of PTnodes to allow breadth-first traversal   */
/**********************************************************************/
typedef struct QueueElement {
  struct PTnode *value;
  struct QueueElement *next;
} QueueElement;
typedef struct Queue {
  struct QueueElement* head;
  struct QueueElement* tail;
} Queue;
/**********************************************************************/

/**********************************************************************/
/* Function prototypes                                                */
/**********************************************************************/
extern void master_process (int, int, char**);
extern void slave_process (int);

extern bool parse_options (int, char**, options_struct*);
extern bool assign_options (char*, char*, options_struct*);
extern void initialize_options (options_struct*);
extern bool validate_options (options_struct*);

extern bool load_initial_coordinates (PTnode*, options_struct*);
extern bool create_root_node (PTnode**, options_struct*);
extern bool initialize_secant (PTnode*, options_struct*);
extern double compute_secant_direction (PTnode*);
extern void assign_predictor_steps (PTnode*, options_struct*);
extern void construct_predictor_nodes (PTnode*, options_struct*);
extern void construct_viable_paths (PTnode*);
extern void choose_viable_paths (PTnode*);
extern void principal_pampac_loop (int, PTnode*, options_struct*);
extern bool write_root_coordinates (PTnode*, options_struct*);

/* These functions deal with the PTnode data structure */
extern PTnode* init_PTnode (int);
extern void delete_tree (PTnode*);
extern void print_PTnode (PTnode*);
extern void print_color (PTnode*, FILE*);
extern void print_tree (PTnode*);
extern void visualize_tree (PTnode*, options_struct*);
extern void assign_processes (PTnode*, int);
extern void assign_depth (PTnode*, int);
extern void assign_color (PTnode*, NodeColors);
extern int count_children (PTnode*);
extern void prune_diverged_nodes (PTnode*, options_struct*);
extern void advance_root_node (PTnode**, options_struct*);

extern void assess_residuals (PTnode*, options_struct*);
extern void compute_corrector_steps (PTnode*, int);
extern void stop_slaves (int);

extern void init_queue (Queue*);
extern void enqueue (Queue*, PTnode*);
extern void dequeue (Queue*);
extern PTnode* front_of_queue (Queue*);
extern int empty_queue (Queue*);

/**********************************************************************/
/* These two functions must be supplied by the user corresponding to  */
/* the desired problem to solve.                                      */
/**********************************************************************/
extern void compute_residual (int, const double*, double*);
extern void single_corrector_step (int, double*, double*);
extern bool write_coordinates (int, double*);

#include "pampac.h"
/**********************************************************************/
/* This routine traverses a tree of PTnodes, generating output in the */
/* dot language (compatible with the GraphViz family of graph tools). */
/**********************************************************************/
void
visualize_tree (PTnode *root, options_struct *opts) {
  int k;
  Queue q;
  PTnode *alpha, *beta;
  char file_name[STRING_LEN];
  FILE *out_file;
  /* Location of the output file determined by parameters in the data
   * structure opts. */
  sprintf (file_name, "%s-%05d.gv",
           opts->tree_base_filename, opts->tree_filename_num);
  out_file = fopen (file_name, "w");
  if (out_file == NULL) {
    printf ("Error Opening %s.\n", file_name);
    exit (1);
  }
  /* Start printing information to dot file */
  fprintf (out_file, "digraph G {\n");
  fprintf (out_file,
           "node [shape=circle,style=filled,fontcolor=black];\n");
  /* Breadth-first traversal of tree */
  init_queue (&q);
  enqueue (&q, root);
  while (!empty_queue (&q)) {
    alpha = front_of_queue (&q);
    dequeue (&q); /* Mark alpha as visited by dequeuing. */

    /* Print node information about alpha into GraphViz file. */
    fprintf (out_file, "%d [label=\"%d\\nnu=%d \", color=",
             alpha->label,alpha->label,alpha->nu);
    print_color (alpha, out_file);
    fprintf (out_file, "];\n");
    /* Append children of alpha just dequeued (if any) to end of queue. */
    if (alpha->child != NULL) {
      for (k=0; k<alpha->max_children; k++) {
        beta = alpha->child[k];
        if (beta != NULL) {
          enqueue (&q, beta);
          /* Print edge information into GraphViz file for the edge
             connecting alpha to beta. */
          fprintf (out_file, "%d->%d [label=\"%-8.4g\"];\n",
                   alpha->label,beta->label,beta->h);
        }
      }
    }
  }
  fprintf (out_file, "}\n");
  fclose (out_file);
  opts->tree_filename_num += 1;
  return;
}

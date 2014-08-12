#include "pampac.h"
/**********************************************************************/
/* Traverse tree (depth-first) assigning colors to all nodes.        */
/**********************************************************************/
void
assign_color (PTnode * node, NodeColors color) {
  int k;
  for (k = 0; k < node->max_children; k++) {
    if (node->child[k] != NULL) {
      assign_color (node->child[k], color);
    }
  }
  node->color = color;
}

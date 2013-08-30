#include "pampac.h"
/**********************************************************************/
/* Given a PTnode, print its color to stdout.                        */
/**********************************************************************/
void
print_color (PTnode *node, FILE *file)
{
  switch (node->color)
    {
    case BLACK:
      fprintf (file, "BLACK");
      break;
    case GREEN:
      fprintf (file, "GREEN");
      break;
    case YELLOW:
      fprintf (file, "YELLOW");
      break;
    case RED:
      fprintf (file, "RED");
      break;
    default:
      break;
    }
}

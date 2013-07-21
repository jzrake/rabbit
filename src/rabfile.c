/* -----------------------------------------------------------------------------
 * FILE:
 *
 * AUTHOR: Jonathan Zrake
 *
 * DESCRIPTION:
 *
 *
 * -----------------------------------------------------------------------------
 */

#define RABBIT_INTERNAL
#include "rabbit.h"


int main(int argc, char **argv)
{
  if (argc == 1) {
    printf("usage: rabfile input.mesh\n");
    return 1;
  }

  rabbit_mesh *mesh = rabbit_mesh_load(argv[1]);
  rabbit_cfg config = mesh->config;
  rabbit_node *node, *tmp_node;
  rabbit_edge *edge, *tmp_edge;
  int n;
  int *I, *V;

  printf("rabbit_cfg:\n");
  printf("  max_depth = %d\n", config.max_depth);
  printf("  doubles_per_node = %d\n", config.doubles_per_node);
  printf("  doubles_per_edge = %d\n", config.doubles_per_edge);

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {
    I = node->index;
    printf("N %d %d %d %d :", I[0], I[1], I[2], I[3]);
    for (n=0; n<mesh->config.doubles_per_node; ++n) {
      printf(" %f", node->data[n]);
    }
    printf("\n");
  }

  HASH_ITER(hh, mesh->edges, edge, tmp_edge) {
    V = edge->vertices;
    printf("E %d %d %d %d %d %d :", V[0], V[1], V[2], V[3], V[4], V[5]);
    for (n=0; n<mesh->config.doubles_per_edge; ++n) {
      printf(" %f", edge->data[n]);
    }
    printf("\n");
  }

  rabbit_mesh_del(mesh);
  return 0;
}

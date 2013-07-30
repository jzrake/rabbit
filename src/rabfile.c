/* -----------------------------------------------------------------------------
 * FILE: rabfile.c
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
  //  rabbit_cfg config = mesh->config;
  rabbit_node *node, *tmp_node;
  //  rabbit_edge *edge, *tmp_edge;
  //  rabbit_geom geom;
  int n;
  //  int *I;
  //  int *V;
  /*
  printf("rabbit_cfg:\n");
  printf("  max_depth = %d\n", config.max_depth);
  printf("  doubles_per_node = %d\n", config.doubles_per_node);
  printf("  doubles_per_face = %d\n", config.doubles_per_face);
  printf("  doubles_per_edge = %d\n", config.doubles_per_edge);

  printf("mesh stats:\n");
  printf("  number of nodes = %d\n", rabbit_mesh_count(mesh, RABBIT_ACTIVE));
  printf("  number of faces = %d\n", rabbit_mesh_count(mesh, RABBIT_FACE));
  printf("  number of edges = %d\n", rabbit_mesh_count(mesh, RABBIT_EDGE));
  */

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {
    //    geom = rabbit_mesh_geom(mesh, node->rnp);
    //    I = geom.index;
    //    printf("N %d %d %d %d :", I[0], I[1], I[2], I[3]);
    printf("%d ", node->rnp[0]);
    for (n=0; n<mesh->config.doubles_per_node; ++n) {
      printf(" %f", node->data[n]);
    }
    printf("\n");
  }

  /*
  HASH_ITER(hh, mesh->edges, edge, tmp_edge) {
    geom = rabbit_mesh_geom(mesh, edge->rnp);
    V = geom.vertices;
    printf("E %d %d %d %d %d %d :", V[0], V[1], V[2], V[3], V[4], V[5]);
    for (n=0; n<mesh->config.doubles_per_edge; ++n) {
      printf(" %f", edge->data[n]);
    }
    printf("\n");
  }
  */

  rabbit_mesh_del(mesh);
  return 0;
}

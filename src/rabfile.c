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

  printf("rabbit_cfg:\n");
  printf("  max_depth = %d\n", config.max_depth);
  printf("  doubles_per_node = %d\n", config.doubles_per_node);
  printf("  doubles_per_edge = %d\n", config.doubles_per_edge);

  rabbit_mesh_del(mesh);

  /*
  while (tpl_unpack(tn, 1) > 0) {
    printf("N %d %d %d %d :", I[0], I[1], I[2], I[3]);
    while (tpl_unpack(tn, 2) > 0) {
      printf(" %+6.4e", node_data_val);
    }
    printf("\n");
  }

  while (tpl_unpack(tn, 3) > 0) {
    printf("E %d %d %d %d %d %d:", V[0], V[1], V[2], V[3], V[4], V[5]);
    while (tpl_unpack(tn, 4) > 0) {
      printf(" %+6.4e", edge_data_val);
    }
    printf("\n");
  }

  tpl_free(tn);
  */

  return 0;
}

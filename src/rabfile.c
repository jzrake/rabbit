#include <stdio.h>
#include "uthash.h"
#include "tpl.h"

typedef struct
{
  int max_depth;
  int doubles_per_node;
  int doubles_per_edge;
} rabbit_cfg;


int main(int argc, char **argv)
{
  if (argc == 1) {
    printf("usage: rabfile input.mesh\n");
    return 1;
  }

  int I[4], V[6];
  rabbit_cfg config_val;
  double node_data_val;
  double edge_data_val;
  tpl_node *tn = tpl_map("S(iii)A(i#A(f))A(i#A(f))",
			 &config_val,     // 0
			 I, 4,            // 1
			 &node_data_val,  // 2
			 V, 6,            // 3
			 &edge_data_val); // 4

  tpl_load(tn, TPL_FILE, argv[1]);
  tpl_unpack(tn, 0);

  printf("rabbit_cfg:\n");
  printf("  max_depth = %d\n", config_val.max_depth);
  printf("  doubles_per_node = %d\n", config_val.doubles_per_node);
  printf("  doubles_per_edge = %d\n", config_val.doubles_per_edge);

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
  return 0;
}

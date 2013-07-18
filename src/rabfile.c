#include <stdio.h>
#include "uthash.h"
#include "tpl.h"

int main(int argc, char **argv)
{
  if (argc == 1) {
    printf("usage: rabfile input.mesh\n");
    return 1;
  }

  int I[4];
  double data_val;
  tpl_node *tn = tpl_map("A(i#A(f))", &I, 4, &data_val);

  tpl_load(tn, TPL_FILE, argv[1]);

  while (tpl_unpack(tn, 1) > 0) {
    printf("%d %d %d %d :", I[0], I[1], I[2], I[3]);
    while (tpl_unpack(tn, 2) > 0) {
      printf(" %+6.4e", data_val);
    }
    printf("\n");
  }

  tpl_free(tn);

  return 0;
}

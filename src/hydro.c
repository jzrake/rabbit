/* -----------------------------------------------------------------------------
 * FILE: hydro.c
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

void build_ghost(rabbit_mesh *mesh)
{
  rabbit_face *face;
  for(face=mesh->faces; face != NULL; face=face->hh.next) {

    int rnpL[3];
    int rnpR[3];
    int h;
    rabbit_node *nL, *nR;
    rabbit_geom geom;

    geom = rabbit_mesh_geom(mesh, face->rnp);
    h = mesh->config.max_depth - geom.index[0] - 1;

    memcpy(rnpL, face->rnp, 3 * sizeof(int));
    memcpy(rnpR, face->rnp, 3 * sizeof(int));

    rnpL[0] -= 1 << h;
    rnpR[0] += 1 << h;

    if (geom.axis == 0) {
      if (rnpL[0] > 0) {
	nL = rabbit_mesh_containing(mesh, rnpL, RABBIT_RNP);
      }
      else {
	nL = rabbit_mesh_getnode(mesh, rnpL, RABBIT_RNP);
      }
      if (rnpR[0] < (1 << mesh->config.max_depth)) {
	nR = rabbit_mesh_containing(mesh, rnpR, RABBIT_RNP);
      }
      else {
	nR = rabbit_mesh_getnode(mesh, rnpR, RABBIT_RNP);
      }

      if (nL == NULL) {
	printf("putting guard at L boundary: %d\n", rnpL[0]);
	rabbit_mesh_putnode(mesh, rnpL, RABBIT_RNP | RABBIT_GHOST);
      }
      if (nR == NULL) {
	printf("putting guard at R boundary: %d\n", rnpR[0]);
	rabbit_mesh_putnode(mesh, rnpR, RABBIT_RNP | RABBIT_GHOST);
      }
    }
  }
  rabbit_mesh_build(mesh);
}

int main(int argc, char **argv)
{
  rabbit_cfg config = { 16, 5, 0, 0 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int d = 3;
  int i;
  int index[4] = { d, 0, 0, 0 };

  for (i=0; i<(1<<d); ++i) {
    index[1] = i;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
  }
  rabbit_mesh_build(mesh);

  for(node=mesh->nodes; node != NULL; node=node->hh.next) {
    double x = ((double) node->rnp[0]) / (1 << mesh->config.max_depth);
    if (x < 0.5) {
      node->data[0] = 1.0;
      node->data[1] = 1.0;
      node->data[2] = 0.0;
      node->data[3] = 0.0;
      node->data[4] = 0.0;
    }
    else {
      node->data[0] = 0.125;
      node->data[1] = 0.1;
      node->data[2] = 0.0;
      node->data[3] = 0.0;
      node->data[4] = 0.0;
    }
  }

  build_ghost(mesh);
  build_ghost(mesh);
  build_ghost(mesh);

  printf("there are %d ghost zones\n", rabbit_mesh_count(mesh, RABBIT_GHOST));

  rabbit_mesh_dump(mesh, "rabbit-1d.mesh");
  rabbit_mesh_del(mesh);
  return 0;
}

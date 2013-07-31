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

#include <math.h>

#define RABBIT_INTERNAL
#include "rabbit.h"
#include "euler1d.h"


void face_get_nodes(rabbit_face *face, rabbit_node **nL, rabbit_node **nR,
		    int put_missing)
{
  int rnpL[3];
  int rnpR[3];
  int h;
  rabbit_geom geom;

  geom = rabbit_mesh_geom(face->mesh, face->rnp);
  h = face->mesh->config.max_depth - geom.index[0] - 1;

  memcpy(rnpL, face->rnp, 3 * sizeof(int));
  memcpy(rnpR, face->rnp, 3 * sizeof(int));

  rnpL[geom.axis] -= 1 << h;
  rnpR[geom.axis] += 1 << h;

  *nL = rabbit_mesh_containing(face->mesh, rnpL, RABBIT_RNP);
  *nR = rabbit_mesh_containing(face->mesh, rnpR, RABBIT_RNP);

  if (put_missing) {
    if (*nL == NULL) {
      *nL = rabbit_mesh_putnode(face->mesh, rnpL, RABBIT_RNP | RABBIT_GHOST);
      memcpy((*nL)->data, (*nR)->data,
	     face->mesh->config.doubles_per_node * sizeof(double));
    }
    if (*nR == NULL) {
      *nR = rabbit_mesh_putnode(face->mesh, rnpR, RABBIT_RNP | RABBIT_GHOST);
      memcpy((*nR)->data, (*nL)->data,
	     face->mesh->config.doubles_per_node * sizeof(double));
    }
  }
}

int face_get_nodes2(rabbit_face *face, rabbit_node **nodes, int num)
{
  int rnp[3];
  int h;
  int n;
  int all_there = 1;
  rabbit_geom geom;

  geom = rabbit_mesh_geom(face->mesh, face->rnp);
  h = face->mesh->config.max_depth - geom.index[0] - 1;

  for (n=0; n<num; ++n) {

    memcpy(rnp, face->rnp, 3 * sizeof(int));
    rnp[geom.axis] -= (2*n + 1) * (1 << h);
    nodes[num - n - 1] = rabbit_mesh_containing(face->mesh, rnp, RABBIT_RNP);
    all_there *= (nodes[num - n - 1] != NULL);

    memcpy(rnp, face->rnp, 3 * sizeof(int));
    rnp[geom.axis] += (2*n + 1) * (1 << h);
    nodes[num + n + 0] = rabbit_mesh_containing(face->mesh, rnp, RABBIT_RNP);
    all_there *= (nodes[num + n + 0] != NULL);

  }

  return all_there;
}

void build_ghost(rabbit_mesh *mesh)
{
  rabbit_face *face;
  rabbit_node *nL, *nR;
  rabbit_geom geom;

  for (face=mesh->faces; face != NULL; face=face->hh.next) {
    geom = rabbit_mesh_geom(mesh, face->rnp);
    if (geom.axis == 0 || geom.axis == 1) {
      face_get_nodes(face, &nL, &nR, 1);
    }
  }
  rabbit_mesh_build(mesh);
}

void onestep(rabbit_mesh *mesh, double dt)
{
  Primitive P[4];
  Primitive Pl, Pr;
  rabbit_face *face, *tmp_face;
  rabbit_geom geomF, geomL, geomR;
  rabbit_node *nodes[4];
  int n, q;

  HASH_ITER(hh, mesh->faces, face, tmp_face) {
    geomF = rabbit_mesh_geom(mesh, face->rnp);
    if (geomF.axis == 0 || geomF.axis == 1) {
      face_get_nodes2(face, nodes, 2);

      if (nodes[1] == NULL) nodes[1] = nodes[2];
      if (nodes[2] == NULL) nodes[2] = nodes[1];
      if (nodes[0] == NULL) nodes[0] = nodes[1];
      if (nodes[3] == NULL) nodes[3] = nodes[2];

      for (n=0; n<4; ++n) {
	P[n] = *((const Primitive*) nodes[n]->data);
      }
      ReconstructStates(&P[1], &Pl, &Pr);
      geomL = rabbit_mesh_geom(mesh, nodes[1]->rnp);
      geomR = rabbit_mesh_geom(mesh, nodes[2]->rnp);
      Conserved F = RiemannSolver(&Pl, &Pr, geomF.axis);
      double dA = 1.0;
      double dVL = 1.0 / (1 << geomL.index[0]);
      double dVR = 1.0 / (1 << geomR.index[0]);

      for (q=0; q<5; ++q) {
	nodes[1]->data[10 + q] -= ((double*)&F)[q] * dt * dA / dVL;
	nodes[2]->data[10 + q] += ((double*)&F)[q] * dt * dA / dVR;
      }
    }
  }
}

void average(rabbit_mesh *mesh, double a, double b)
{
  rabbit_node *node, *tmp_node;
  int q;

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {
    if ((node->flags & RABBIT_GHOST) == 0) {
      double *U0 = &node->data[5];
      double *U1 = &node->data[10];
      for (q=0; q<5; ++q) {
	U1[q] = a * U0[q] + b * U1[q];
      }
      const Primitive P = ConsToPrim((const Conserved*) U1);
      memcpy(&node->data[0], &P, 5 * sizeof(double));
    }
  }
}

void advance(rabbit_mesh *mesh, double dt)
{
  rabbit_node *node, *tmp_node;

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {
    const Primitive *P = (const Primitive*) &node->data[0];
    const Conserved U = PrimToCons(P);
    memcpy(&node->data[ 5], &U, 5 * sizeof(double));
    memcpy(&node->data[10], &U, 5 * sizeof(double));
  }

  onestep(mesh, dt);
  average(mesh, 0.0, 1.0);
}

void advanceRK3(rabbit_mesh *mesh, double dt)
{
  rabbit_node *node, *tmp_node;

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {
    const Primitive *P = (const Primitive*) &node->data[0];
    const Conserved U = PrimToCons(P);
    memcpy(&node->data[ 5], &U, 5 * sizeof(double));
    memcpy(&node->data[10], &U, 5 * sizeof(double));
  }

  onestep(mesh, dt);
  average(mesh, 0./1, 1./1);

  onestep(mesh, dt);
  average(mesh, 3./4, 1./4);

  onestep(mesh, dt);
  average(mesh, 1./3, 2./3);
}

int main(int argc, char **argv)
{
  rabbit_cfg config = { 16, 15, 0, 0 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int d = 7;
  int i, j;
  int index[4] = { d, 0, 0, 0 };

  if (1) {
    for (i=2; i<(1<<d)-2; ++i) {
      for (j=2; j<(1<<d)-2; ++j) {
	index[0] = d;
	index[1] = i;
	index[2] = j;
	rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
      }
    }
  }

  rabbit_mesh_build(mesh);

  for (node=mesh->nodes; node != NULL; node=node->hh.next) {
    double x = ((double) node->rnp[0]) / (1 << mesh->config.max_depth) - 0.5;
    double y = ((double) node->rnp[1]) / (1 << mesh->config.max_depth) - 0.5;
    double r = sqrt(x*x + y*y);

    if (r < 0.25) {
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

  //  build_ghost(mesh);
  //  build_ghost(mesh);
  printf("there are %d ghost zones\n", rabbit_mesh_count(mesh, RABBIT_GHOST));


  double t = 0.0;
  double dt = 0.00075;

  while (t < 0.025) {
    advanceRK3(mesh, dt);
    t += dt;
    printf("t=%4.3f\n", t);
  }


  rabbit_mesh_dump(mesh, "rabbit-2d.mesh");
  rabbit_mesh_del(mesh);
  return 0;
}

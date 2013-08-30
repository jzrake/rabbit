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


/* Maximum value for a rational number position */
#define max_rnp(M) (1 << (M)->config.max_depth)

static const double PlmTheta = 2.0;
static double minval3(double a, double b, double c)
{
  return (a<b) ? ((a<c) ? a : c) : ((b<c) ? b : c);
}
static double sign(double x)
{
  return (x > 0.0) - (x < 0.0);
}


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
      printf("putting guard at L boundary: %d\n", (*nL)->rnp[0]);
      memcpy((*nL)->data, (*nR)->data,
	     face->mesh->config.doubles_per_node * sizeof(double));
    }
    if (*nR == NULL) {
      *nR = rabbit_mesh_putnode(face->mesh, rnpR, RABBIT_RNP | RABBIT_GHOST);
      printf("putting guard at R boundary: %d\n", (*nR)->rnp[0]);
      memcpy((*nR)->data, (*nL)->data,
	     face->mesh->config.doubles_per_node * sizeof(double));
    }
  }
}

void build_ghost(rabbit_mesh *mesh)
{
  rabbit_face *face;
  rabbit_node *nL, *nR;
  rabbit_geom geom;

  for (face=mesh->faces; face != NULL; face=face->hh.next) {
    geom = rabbit_mesh_geom(mesh, face->rnp);
    if (geom.axis == 0) {
      face_get_nodes(face, &nL, &nR, 1);
    }
  }
  rabbit_mesh_build(mesh);
}

void evalgrad(rabbit_mesh *mesh)
{
  rabbit_node *node, *tmp_node;
  rabbit_node *nodeL;
  rabbit_node *nodeR;
  rabbit_geom geomL;
  rabbit_geom geomR;
  int q;

  rabbit_mesh_sort(mesh);

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {

    geomL = geomR = rabbit_mesh_geom(mesh, node->rnp);
    geomL.index[1] -= 1;
    geomR.index[1] += 1;

    nodeL = rabbit_mesh_containing(mesh, geomL.index, RABBIT_INDEX);
    nodeR = rabbit_mesh_containing(mesh, geomR.index, RABBIT_INDEX);

    if (nodeL == NULL) {
      int size;
      nodeL = rabbit_mesh_contains(mesh, geomL.index, RABBIT_INDEX, &size);
      if (nodeL) {
	nodeL = nodeL->hh.next; // need the right node in the subtree
      }
    }
    if (nodeR == NULL) {
      int size;
      nodeR = rabbit_mesh_contains(mesh, geomR.index, RABBIT_INDEX, &size);
    }

    if (!(nodeL && nodeR)) {
      for (q=0; q<5; ++q) {
	node->data[15 + q] = 0.0;
      }
      node->flags |= RABBIT_GHOST;
    }
    else {
      double xL = (double) nodeL->rnp[0] / max_rnp(mesh);
      double x0 = (double) node ->rnp[0] / max_rnp(mesh);
      double xR = (double) nodeR->rnp[0] / max_rnp(mesh);

      for (q=0; q<5; ++q) {
	double yL = (double) nodeL->data[q];
	double y0 = (double) node ->data[q];
	double yR = (double) nodeR->data[q];
	
	double a = (yL - y0) / (xL - x0) * PlmTheta;
	double b = (yL - yR) / (xL - xR);
	double c = (y0 - yR) / (x0 - xR) * PlmTheta;
	
	double sa = sign(a), sb = sign(b), sc = sign(c);
	double fa = fabs(a), fb = fabs(b), fc = fabs(c);
	double g = 0.25 * fabs(sa + sb) * (sa + sc) * minval3(fa, fb, fc);
	
	node->data[15 + q] = g;
      }
    }
  }
}

void onestep(rabbit_mesh *mesh, double dt)
{
  Primitive Pl, Pr;
  rabbit_face *face, *tmp_face;
  rabbit_geom geomF, geomL, geomR;
  rabbit_node *nodeL, *nodeR;
  double gL, gR;
  double dL, dR;
  int q;

  HASH_ITER(hh, mesh->faces, face, tmp_face) {

    geomF = rabbit_mesh_geom(mesh, face->rnp);

    if (geomF.axis == 0) {

      face_get_nodes(face, &nodeL, &nodeR, 0);

      if (nodeL == NULL || nodeR == NULL) {
	//	printf("on boundary\n");
	continue;
      }

      geomL = rabbit_mesh_geom(mesh, nodeL->rnp);
      geomR = rabbit_mesh_geom(mesh, nodeR->rnp);

      for (q=0; q<5; ++q) {

	dL = (double) (face->rnp[0] - nodeL->rnp[0]) / max_rnp(mesh);
	dR = (double) (face->rnp[0] - nodeR->rnp[0]) / max_rnp(mesh);
	gL = nodeL->data[15 + q];
	gR = nodeR->data[15 + q];

	((double*) &Pl)[q] = nodeL->data[q] + dL * gL;
	((double*) &Pr)[q] = nodeR->data[q] + dR * gR;
      }

      Conserved F = RiemannSolver(&Pl, &Pr, 0);

      double dA = 1.0;
      double dVL = 1.0 / (1 << geomL.index[0]);
      double dVR = 1.0 / (1 << geomR.index[0]);

      for (q=0; q<5; ++q) {
	nodeL->data[10 + q] -= ((double*)&F)[q] * dt * dA / dVL;
	nodeR->data[10 + q] += ((double*)&F)[q] * dt * dA / dVR;
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

void prepare_timestep(rabbit_mesh *mesh)
{
  rabbit_node *node, *tmp_node;

  HASH_ITER(hh, mesh->nodes, node, tmp_node) {

    if (node->flags & RABBIT_GHOST) continue;

    const Primitive *P = (const Primitive*) &node->data[0];
    const Conserved U = PrimToCons(P);

    memcpy(&node->data[ 5], &U, 5 * sizeof(double));
    memcpy(&node->data[10], &U, 5 * sizeof(double));
  }
}

void advance(rabbit_mesh *mesh, double dt)
{
  prepare_timestep(mesh);

  onestep(mesh, dt);
  average(mesh, 0.0, 1.0);
  evalgrad(mesh);
}

void advanceRK3(rabbit_mesh *mesh, double dt)
{
  prepare_timestep(mesh);

  evalgrad(mesh);
  onestep(mesh, dt);
  average(mesh, 0./1, 1./1);

  evalgrad(mesh);
  onestep(mesh, dt);
  average(mesh, 3./4, 1./4);

  evalgrad(mesh);
  onestep(mesh, dt);
  average(mesh, 1./3, 2./3);
}

int main(int argc, char **argv)
{
  rabbit_cfg config = { 11, 20, 0, 0 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int d = 9;
  int i;
  int index[4] = { d, 0, 0, 0 };

  if (0) {
    for (i=8; i<(1<<d)-8; ++i) {
      index[1] = i;
      rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    }
  }
  else if (1) {
    for (i=0; i<2*(1<<d)/8; ++i) {
      index[0] = d - 1;
      index[1] = i;
      rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    }
    for (i=2*(1<<d)/4; i<(1<<d); ++i) {
      index[0] = d;
      index[1] = i;
      rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    }
  }
  rabbit_mesh_build(mesh);

  for (node=mesh->nodes; node != NULL; node=node->hh.next) {
    double x = ((double) node->rnp[0]) / max_rnp(mesh);
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

  double t = 0.0;
  double dt = 0.00075;

  while (t < 0.20) {
    advanceRK3(mesh, dt);
    t += dt;
    printf("t=%4.3f\n", t);
  }

  rabbit_mesh_dump(mesh, "rabbit-1d.mesh");
  rabbit_mesh_del(mesh);
  return 0;
}

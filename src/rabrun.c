#include <stdio.h>
#include "uthash.h"
#include "tpl.h"

#define DEBUG 1
#define DEBUG_MSG(m) if (DEBUG) { printf("%s %s\n", __FUNCTION__, m); }

#define RABBIT_ANY      (1 << 0)
#define RABBIT_ACTIVE   (1 << 1)
#define RABBIT_GHOST    (1 << 2)
#define RABBIT_FORCE    (1 << 3)



typedef struct rabbit_mesh rabbit_mesh;
typedef struct rabbit_node rabbit_node;
typedef struct rabbit_face rabbit_face;
typedef struct rabbit_edge rabbit_edge;
typedef struct
{
  int rank;
  int max_depth;
  int doubles_per_node;
  int doubles_per_edge;
} rabbit_cfg;


void         rabbit_mesh_del(rabbit_mesh *M);
rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg);
rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int index[4], int flags);
rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_containing(rabbit_mesh *M, int index[4]);
int          rabbit_mesh_count(rabbit_mesh *M, int flags);
void         rabbit_mesh_build(rabbit_mesh *M);
void         rabbit_mesh_dump(rabbit_mesh *M, char *fname);


struct rabbit_mesh {
  rabbit_cfg config;
  rabbit_node *nodes;
  rabbit_edge *edges;
} ;

struct rabbit_node {
  int index[4]; // (depth, i, j, k)
  int flags;
  double *data;
  rabbit_mesh *mesh;
  UT_hash_handle hh;
} ;

struct rabbit_face {
  rabbit_node nodes[2];
} ;

struct rabbit_edge {
  int vertices[6];
  double *data;
  UT_hash_handle hh;
} ;


rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg)
{
  rabbit_mesh *M = (rabbit_mesh*) malloc(sizeof(rabbit_mesh));
  M->config = cfg;
  M->nodes = NULL;
  M->edges = NULL;
  return M;
}

void rabbit_mesh_del(rabbit_mesh *M)
{
  rabbit_node *node, *tmp_node;
  rabbit_edge *edge, *tmp_edge;

  HASH_ITER(hh, M->nodes, node, tmp_node) {
    rabbit_mesh_delnode(M, node->index);
  }
  HASH_ITER(hh, M->edges, edge, tmp_edge) {

    HASH_DEL(M->edges, edge);

    printf("removing edge, there are now %d edges\n",
           HASH_CNT(hh, M->edges));

    free(edge->data);
    free(edge);
  }
  free(M);
}

rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int index[4], int flags)
{
  rabbit_node *node;

  HASH_FIND(hh, M->nodes, index, 4 * sizeof(int), node);

  if (node != NULL) {
    return node;
  }
  else {
    node = (rabbit_node*) malloc(sizeof(rabbit_node));
    node->flags = flags & (RABBIT_ACTIVE | RABBIT_GHOST);
    node->mesh = M;
    node->data = (double*) calloc(M->config.doubles_per_node, sizeof(double));
    memcpy(node->index, index, 4 * sizeof(int));
    HASH_ADD(hh, M->nodes, index, 4 * sizeof(int), node);
    return node;
  }
}

rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int index[4])
{
  rabbit_node *node;

  HASH_FIND(hh, M->nodes, index, 4 * sizeof(int), node);

  return node;
}

rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int index[4])
{
  rabbit_node *node;

  HASH_FIND(hh, M->nodes, index, 4 * sizeof(int), node);
  HASH_DEL(M->nodes, node);
  free(node->data);
  free(node);

  return NULL;
}

rabbit_node *rabbit_mesh_containing(rabbit_mesh *M, int index[4])
{
  rabbit_node *node;
  int I[4];

  if (index[0] < 0) {
    return NULL;
  }
  else if ((node = rabbit_mesh_getnode(M, index))) {
    return node;
  }
  else {
    I[0] = index[0] - 1;
    I[1] = index[1] / 2;
    I[2] = index[2] / 2;
    I[3] = index[3] / 2;
    return rabbit_mesh_containing(M, I);
  }
}

int rabbit_mesh_count(rabbit_mesh *M, int flags)
{
  rabbit_node *node, *tmp;
  int count = 0;

  if (flags & RABBIT_ANY) {
    return HASH_CNT(hh, M->nodes);
  }
  else {
    HASH_ITER(hh, M->nodes, node, tmp) {
      if (node->flags & flags) {
        count += 1;
      }
    }
    return count;
  }
}

void rabbit_mesh_build(rabbit_mesh *M)
{
  rabbit_node *node, *tmp;
  rabbit_edge *edge = NULL;
  rabbit_edge *existing_edge = NULL;

  int vertices[8][3];
  int start[3][4] = {{0, 2, 4, 6},
                     {0, 4, 1, 5},
                     {0, 1, 2, 3}};
  int jumps[3] = { 1, 2, 4 };

  /* ------------------
   * X      Y      Z
   * ------------------
   * 0->1   0->2   0->4
   * 2->3   4->6   1->5
   * 4->5   1->3   2->6
   * 6->7   5->7   3->7
   * ------------------
   */


  HASH_ITER(hh, M->nodes, node, tmp) {


    int n; // starting node counter, [0, 3)
    int a, ai; // axis counter, [0,4)
    int d = node->index[0]; // depth
    int f = M->config.max_depth - d + 1; // rational index log2 denominator


    /* set the index (i,j,k) of each of this node's 8 vertices */

    for (n=0; n<8; ++n) {

      vertices[n][0] = (node->index[1] + (n >> 0)) << f;
      vertices[n][1] = (node->index[2] + (n >> 1)) << f;
      vertices[n][2] = (node->index[3] + (n >> 2)) << f;

    }


    /* loop over 3 axes 'a' */

    for (a=0; a<3; ++a) {


      /* loop over 4 starting vertices 's' */

      for (n=0; n<4; ++n) {

        HASH_FIND(hh, M->edges, vertices, 6 * sizeof(int), existing_edge);

        if (existing_edge == NULL) {

	  int v0 = start[a][n];
	  int v1 = start[a][n] + jumps[a];

          edge = (rabbit_edge*) malloc(sizeof(rabbit_edge));
          edge->data = (double*) calloc(M->config.doubles_per_edge, sizeof(double));

	  for (ai=0; ai<3; ++ai) {
	    edge->vertices[ai+0] = vertices[v0][a];
	    edge->vertices[ai+3] = vertices[v1][a];
	  }

	  HASH_ADD(hh, M->edges, vertices, 6 * sizeof(int), edge);

          printf("there are now %d edges\n", HASH_CNT(hh, M->edges));

        }
      }
    }
  }
}

void rabbit_mesh_dump(rabbit_mesh *M, char *fname)
{
  int n, I[4];
  double data_val;
  rabbit_node *node, *tmp;
  tpl_node *tn = tpl_map("A(i#A(f))", &I, 4, &data_val);

  HASH_ITER(hh, M->nodes, node, tmp) {

    memcpy(I, node->index, 4 * sizeof(int));
    tpl_pack(tn, 1);

    for (n=0; n<M->config.doubles_per_node; ++n) {
      data_val = node->data[n];
      tpl_pack(tn, 2);
    }
  }

  tpl_dump(tn, TPL_FILE, fname);
  tpl_free(tn);
}


int main()
{
  rabbit_cfg config = { 3, 10, 4, 4 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int I[4] = { 0, 0, 0, 0 };
  int i,j,k;
  int D = 0;

  for (i=0; i<(1<<D); ++i) {
    for (j=0; j<(1<<D); ++j) {
      for (k=0; k<(1<<D); ++k) {
        I[0] = D;
        I[1] = i;
        I[2] = j;
        I[3] = k;
        node = rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
      }
    }
  }

  rabbit_mesh_build(mesh);

  printf("there are %d total nodes\n", rabbit_mesh_count(mesh, RABBIT_ANY));

  rabbit_mesh_dump(mesh, "rabbit.mesh");

  rabbit_mesh_del(mesh);

  return 0;
}

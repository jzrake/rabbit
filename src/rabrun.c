#include <stdio.h>
#include "uthash.h"
#include "tpl.h"

#define RABBIT_NODE     (1 << 0)
#define RABBIT_GHOST    (1 << 1)
#define RABBIT_FORCE    (1 << 2)
#define RABBIT_CONTAIN  (1 << 3)

typedef struct rabbit_mesh rabbit_mesh;
typedef struct rabbit_node rabbit_node;
typedef struct rabbit_face rabbit_face;
typedef struct rabbit_edge rabbit_edge;
typedef struct
{
  int rank;
  int max_depth;
  int doubles_per_node;
} rabbit_cfg;


rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg);
rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int index[4]);
void rabbit_mesh_dump(rabbit_mesh *M, char *fname);

struct rabbit_mesh {
  rabbit_cfg config;
  rabbit_node *nodes;
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
  int vertex0[3];
  int vertex1[3];
} ;


rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg)
{
  rabbit_mesh *M = (rabbit_mesh*) malloc(sizeof(rabbit_mesh));
  M->config = cfg;
  M->nodes = NULL;
  return M;
}

void rabbit_mesh_del(rabbit_mesh *M)
{
  struct rabbit_node *node, *tmp;
  HASH_ITER(hh, M->nodes, node, tmp) {
    rabbit_mesh_delnode(M, node->index);
  }
  free(M);
}

rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int index[4])
{
  rabbit_node *node;
  HASH_FIND(hh, M->nodes, index, 4 * sizeof(int), node);
  if (node != NULL) {
    return node;
  }
  else {
    node = (rabbit_node*) malloc(sizeof(rabbit_node));
    node->flags = RABBIT_NODE;
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

int rabbit_mesh_count(rabbit_mesh *M)
{
  return HASH_CNT(hh, M->nodes);
}

void rabbit_mesh_serialize(rabbit_mesh *M, char *fname)
{
  int n, I[4];
  double data_val;
  tpl_node *tn = tpl_map("A(i#A(f))", &I, 4, &data_val);
  struct rabbit_node *node, *tmp;
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
  rabbit_cfg config = { 3, 10, 4 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  int I[4] = { 0, 0, 0, 0 };
  int i,j;
  int D = 6;

  for (i=0; i<(1<<D); ++i) {
    for (j=0; j<(1<<D); ++j) {
      I[0] = D;
      I[1] = i;
      I[2] = j;
      rabbit_node *N = rabbit_mesh_putnode(mesh, I);
      N->data[0] = 1.0;
    }
  }

  printf("there are %d total nodes\n", rabbit_mesh_count(mesh));
  rabbit_mesh_serialize(mesh, "rabbit.mesh");
  rabbit_mesh_del(mesh);
  return 0;
}

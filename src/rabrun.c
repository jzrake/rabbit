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

#include <stdio.h>

/* debug level options */
#define ALWAYS 3
#define SOMETIMES 2
#define ALMOST_NEVER 1
#define NEVER 0

/* print messages how often? */
#define PRINT_MESSAGES ALMOST_NEVER

#define MSG(level, format, ...) do {            \
    if (level < PRINT_MESSAGES) {               \
      fprintf(stderr, "[%s]$ ",  __FUNCTION__); \
      fprintf(stderr, format, __VA_ARGS__);     \
      fprintf(stderr, "\n");                    \
    }                                           \
  } while (0)                                   \

#include <time.h>
#define TIME(cmd) do {                                  \
    clock_t start = clock();                            \
    cmd;                                                \
    printf("[%s]$ %s took %5.4f ms\n", __FUNCTION__,    \
           #cmd, 1e3*(clock() - start)/CLOCKS_PER_SEC); \
  } while (0)


/* assertion macro */
#include <assert.h>
#define NOHUP 0 // continue even if an assertion fails
#define ASSERT_MSG "[assertion:%s]$ %s == %d : %d\n"
#define ASSERTEQ(E,v)printf(ASSERT_MSG,__FUNCTION__,#E,v,E);assert(E==v||NOHUP);


/* ----------------------------------
 * Hash functions available in uthash
 * ----------------------------------
 * JEN Jenkins (default)
 * BER Bernstein
 * SAX Shift-Add-Xor
 * OAT One-at-a-time
 * FNV Fowler/Noll/Vo
 * SFH Paul Hsieh
 * MUR MurmurHash v3 (see note)
 */
#define HASH_FUNCTION HASH_SAX
#include "uthash.h"

/* tpl file format header */
#include "tpl.h"


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


/*
 *
 * private functions
 *
 */

static uint64_t node_preorder_label(rabbit_node *node);
static uint64_t interleave_bits3(uint64_t a, uint64_t b, uint64_t c);
static int      edge_contiguous_compare(rabbit_edge *A, rabbit_edge *B);
static int      edge_contains(rabbit_edge *A, rabbit_edge *B);

/*
 * Return the size of a tree of max depth depth n and branching ratio m=2^r
 * tree_size_atlevel = (m^(n+1) - 1) / (m - 1)
 */
#define tree_size_atlevel(r, n) ((1 << ((r)*((n)+1))) - 1) / ((1 << (r)) - 1)




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
    MSG(2, "removing edge %d", HASH_CNT(hh, M->edges));

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
    MSG(1, "added node with preorder label %"PRIu64,
        node_preorder_label(node));

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
  rabbit_node *node, *tmp_node;
  rabbit_edge *edge = NULL, *tmp_edge, *last_edge;
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


  HASH_ITER(hh, M->nodes, node, tmp_node) {

    int n; // starting node counter, [0, 4)
    int a, ai; // axis counter, [0, 3)
    int d = node->index[0]; // depth
    int f = M->config.max_depth - d + 2; // rational index log2 denominator

    /* set the index (i,j,k) of each of this node's 8 vertices */

    vertices[0][0] = (node->index[1] + 0) << f;
    vertices[0][1] = (node->index[2] + 0) << f;
    vertices[0][2] = (node->index[3] + 0) << f;
    vertices[1][0] = (node->index[1] + 0) << f;
    vertices[1][1] = (node->index[2] + 0) << f;
    vertices[1][2] = (node->index[3] + 1) << f;
    vertices[2][0] = (node->index[1] + 0) << f;
    vertices[2][1] = (node->index[2] + 1) << f;
    vertices[2][2] = (node->index[3] + 0) << f;
    vertices[3][0] = (node->index[1] + 0) << f;
    vertices[3][1] = (node->index[2] + 1) << f;
    vertices[3][2] = (node->index[3] + 1) << f;
    vertices[4][0] = (node->index[1] + 1) << f;
    vertices[4][1] = (node->index[2] + 0) << f;
    vertices[4][2] = (node->index[3] + 0) << f;
    vertices[5][0] = (node->index[1] + 1) << f;
    vertices[5][1] = (node->index[2] + 0) << f;
    vertices[5][2] = (node->index[3] + 1) << f;
    vertices[6][0] = (node->index[1] + 1) << f;
    vertices[6][1] = (node->index[2] + 1) << f;
    vertices[6][2] = (node->index[3] + 0) << f;
    vertices[7][0] = (node->index[1] + 1) << f;
    vertices[7][1] = (node->index[2] + 1) << f;
    vertices[7][2] = (node->index[3] + 1) << f;


    /* loop over 3 axes 'a' */

    for (a=0; a<3; ++a) {


      /* loop over 4 starting vertices 's' */

      for (n=0; n<4; ++n) {

        int v0 = start[a][n];
        int v1 = start[a][n] + jumps[a];
        int edge_vertices[6];

        for (ai=0; ai<3; ++ai) {
          edge_vertices[ai+0] = vertices[v0][ai];
          edge_vertices[ai+3] = vertices[v1][ai];
        }

        HASH_FIND(hh, M->edges, edge_vertices, 6 * sizeof(int), existing_edge);


        if (existing_edge != NULL) {

          MSG(2, "edge exists : [%d %d %d] -> [%d %d %d]",
              edge_vertices[0], edge_vertices[1], edge_vertices[2],
              edge_vertices[3], edge_vertices[4], edge_vertices[5]);
        }
        else {

          edge = (rabbit_edge*) malloc(sizeof(rabbit_edge));
          edge->data = (double*) calloc(M->config.doubles_per_edge,
                                        sizeof(double));

          memcpy(edge->vertices, edge_vertices, 6 * sizeof(int));

          HASH_ADD(hh, M->edges, vertices, 6 * sizeof(int), edge);

          MSG(2, "adding edge %d [%d %d %d] -> [%d %d %d] (v%d -> v%d)",
              HASH_CNT(hh, M->edges),
              edge->vertices[0], edge->vertices[1], edge->vertices[2],
              edge->vertices[3], edge->vertices[4], edge->vertices[5],
	      v0, v1);
        }
      }
    }
  }

  HASH_SRT(hh, M->edges, edge_contiguous_compare);

  last_edge = NULL;
  edge = NULL;

  int iter = 0;
  int removed = 0;


  HASH_ITER(hh, M->edges, edge, tmp_edge) {

    MSG(2, "checking edge %d", iter++);

    if (last_edge != NULL) {

      if (edge_contains(last_edge, edge)) {

        HASH_DEL(M->edges, last_edge);

        MSG(2, "removing duplicate edge %d", removed++);

        free(edge->data);
        free(edge);
      }
    }
    last_edge = edge;
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

int edge_contiguous_compare(rabbit_edge *A, rabbit_edge *B)
{
  int n, len_a, len_b, axis_a, axis_b, depth_a=0, depth_b=0;
  int ax0, ax1, ax2;

  for (n=0; n<3; ++n) {

    len_a = A->vertices[n+3] - A->vertices[n];
    len_b = B->vertices[n+3] - B->vertices[n];

    if (len_a != 0) {
      axis_a = n;
      while (len_a >>= 1) ++depth_a;
    }

    if (len_b != 0) {
      axis_b = n;
      while (len_b >>= 1) ++depth_b;
    }
  }

  ax0 = axis_a;
  ax1 = (ax0 + 1) % 3; // only used if axis_a == axis_b
  ax2 = (ax0 + 2) % 3;

  /* orientation (x, y, z) - directed */
  if (axis_a != axis_b) {
    return axis_b - axis_a;
  }

  /* coordinate of next axis */
  if (A->vertices[ax1] != B->vertices[ax1]) {
    return B->vertices[ax1] - A->vertices[ax1];
  }

  /* coordinate of next axis */
  if (A->vertices[ax2] != B->vertices[ax2]) {
    return B->vertices[ax2] - A->vertices[ax2];
  }

  /* left endpoint of segment along its own axis */
  if (A->vertices[ax0] != B->vertices[ax0]) {
    return B->vertices[ax0] - A->vertices[ax0];
  }

  /* depth of segment */
  if (depth_a != depth_b) {
    return depth_b - depth_a;
  }

  return 0;
}

int edge_contains(rabbit_edge *A, rabbit_edge *B)
/*
 * Return true if edge A contains edge B
 */
{
  int n, len_a, len_b, axis_a, axis_b;
  int ax0, ax1, ax2;

  for (n=0; n<3; ++n) {

    len_a = A->vertices[n+3] - A->vertices[n];
    len_b = B->vertices[n+3] - B->vertices[n];

    if (len_a != 0) {
      axis_a = n;
    }

    if (len_b != 0) {
      axis_b = n;
    }
  }

  ax0 = axis_a;
  ax1 = (ax0 + 1) % 3; // only used if axis_a == axis_b
  ax2 = (ax0 + 2) % 3;

  /* different orientation? */
  if (axis_a != axis_b) {
    return 0;
  }

  /* not co-linear? */
  if (A->vertices[ax1] != B->vertices[ax1]) {
    return 0;
  }
  if (A->vertices[ax2] != B->vertices[ax2]) {
    return 0;
  }

  return (A->vertices[ax0+0] <= B->vertices[ax0+0] &&
          A->vertices[ax0+3] >= B->vertices[ax0+3]);
}

int node_preorder_compare(rabbit_node *a, rabbit_node *b)
{
  return node_preorder_label(a) - node_preorder_label(b);
}

uint64_t node_preorder_label(rabbit_node *node)
/*
 * Return the order in which a given node is visited in a preorder traversal of
 * a fully fleshed out tree having arbitrary max_depth and branching ratio m=8
 */
{
  int r = 3;
  int m = 1 << r; // branching ratio
  int d, n, h, nb, sd, Md, adding, label=0;
  uint64_t index = interleave_bits3(node->index[1],
                                    node->index[2],
                                    node->index[3]);

  for (d=0; d<node->index[0]; ++d) {

    n = node->index[0] - d - 1; // bit
    h = node->mesh->config.max_depth - d - 1;
    nb = 1 << (3*n);
    sd = (index / nb) % m;
    Md = tree_size_atlevel(r, h);
    adding = sd * Md + 1;
    label += adding;

  }

  return label;
}


uint64_t interleave_bits3(uint64_t a, uint64_t b, uint64_t c)
/*
 * Create a 64-bit integer by interleaving the bits of the 21-bit integers a, b,
 * and c
 */
{
  int n;
  uint64_t label = 0;

  for (n=0; n<21; ++n) {
    label += ((a >> n) & 1) << (3*n + 0);
    label += ((b >> n) & 1) << (3*n + 1);
    label += ((c >> n) & 1) << (3*n + 2);
  }

  return label;
}

static void sanity_tests()
{
  /* 1d trees, m=2 */
  ASSERTEQ(tree_size_atlevel(1, 0), 1);
  ASSERTEQ(tree_size_atlevel(1, 1), 3);
  ASSERTEQ(tree_size_atlevel(1, 2), 7);

  /* 2d trees, m=4 */
  ASSERTEQ(tree_size_atlevel(2, 0), 1);
  ASSERTEQ(tree_size_atlevel(2, 1), 5);
  ASSERTEQ(tree_size_atlevel(2, 2), 21);

  /* 3d trees, m=8 */
  ASSERTEQ(tree_size_atlevel(3, 0), 1);
  ASSERTEQ(tree_size_atlevel(3, 1), 9);
  ASSERTEQ(tree_size_atlevel(3, 2), 73);
}

int main()
{
  sanity_tests();

  rabbit_cfg config = { 10, 4, 4 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int I[4] = { 0, 0, 0, 0 };
  int i,j,k;
  int D = 5;

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

  TIME(
       rabbit_mesh_build(mesh)
       );

  TIME(
       HASH_SRT(hh, mesh->nodes, node_preorder_compare)
       );

  TIME(
       HASH_SRT(hh, mesh->edges, edge_contiguous_compare)
       );

  MSG(0, "there are %d total nodes",
      rabbit_mesh_count(mesh, RABBIT_ANY));

  rabbit_mesh_dump(mesh, "rabbit.mesh");

  rabbit_mesh_del(mesh);

  return 0;
}

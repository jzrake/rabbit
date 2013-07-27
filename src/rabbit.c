/* -----------------------------------------------------------------------------
 * FILE: rabbit.c
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

/*
 *
 * private functions
 *
 */
static uint64_t preorder_label(int index[4], int max_depth, int r);
static uint64_t node_preorder_label(rabbit_node *node);
static uint64_t interleave_bits2(uint64_t a, uint64_t b);
static uint64_t interleave_bits3(uint64_t a, uint64_t b, uint64_t c);
//static int      node_preorder_compare(rabbit_node *A, rabbit_node *B);
static int      face_contiguous_compare(rabbit_face *A, rabbit_face *B);
static int      edge_contiguous_compare(rabbit_edge *A, rabbit_edge *B);
static int      face_contains(rabbit_face *A, rabbit_face *B);
static int      edge_contains(rabbit_edge *A, rabbit_edge *B);

/*
 * Return the size of a tree of max depth depth n and branching ratio m=2^r
 * tree_size_atlevel = (m^(n+1) - 1) / (m - 1)
 */
#define tree_size_atlevel(r, n) ((1 << ((r)*((n)+1))) - 1) / ((1 << (r)) - 1)


rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg)
{
  rabbit_mesh *M = (rabbit_mesh*) malloc(sizeof(rabbit_mesh));
  M->config = cfg;
  M->nodes = NULL;
  M->faces = NULL;
  M->edges = NULL;
  return M;
}

void rabbit_mesh_del(rabbit_mesh *M)
{
  rabbit_node *node, *tmp_node;
  rabbit_face *face, *tmp_face;
  rabbit_edge *edge, *tmp_edge;

  HASH_ITER(hh, M->nodes, node, tmp_node) {
    HASH_DEL(M->nodes, node);
    free(node->data);
    free(node);
  }

  HASH_ITER(hh, M->faces, face, tmp_face) {
    HASH_DEL(M->faces, face);
    free(face->data);
    free(face);
  }

  HASH_ITER(hh, M->edges, edge, tmp_edge) {
    HASH_DEL(M->edges, edge);
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
    MSG(1, "added node with preorder label %"PRIu64, node_preorder_label(node));

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

  if (node != NULL) {
    HASH_DEL(M->nodes, node);
    free(node->data);
    free(node);
  }

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

rabbit_geom rabbit_mesh_geom(rabbit_mesh *M, int rnp[3])
{
  rabbit_geom geom;
  int D = M->config.max_depth;
  int H[3] = { 0, 0, 0 }; // height
  int i, j, k, m;
  int ax0, ax1, ax2;

  /* find the position of the least significant active bit in the rational
     number position, the "height" */
  while (((rnp[0] >> H[0]) & 1) == 0 && H[0] < D) ++H[0];
  while (((rnp[1] >> H[1]) & 1) == 0 && H[1] < D) ++H[1];
  while (((rnp[2] >> H[2]) & 1) == 0 && H[2] < D) ++H[2];

  /* if the height is the same along every axis then this is a node */
  if (H[0] == H[1] && H[1] == H[2]) {
    geom.type = RABBIT_NODE;
    geom.axis = -1; /* no orientation */
    geom.index[0] = D - H[0] - 1;
    geom.index[1] = ((rnp[0] >> H[0]) - 1) >> 1;
    geom.index[2] = ((rnp[1] >> H[0]) - 1) >> 1;
    geom.index[3] = ((rnp[2] >> H[0]) - 1) >> 1;
  }

  /* otherwise it's a face or edge */
  else {

    if (H[1] == H[2]) { ax0 = 0; ax1 = 1; ax2 = 2; }
    if (H[2] == H[0]) { ax0 = 1; ax1 = 2; ax2 = 0; }
    if (H[0] == H[1]) { ax0 = 2; ax1 = 0; ax2 = 1; }

    /* otherwise this is an edge or a face - for an edge the double height value
       is larger than the remaining height value */

    if (H[ax0] > H[ax1]) {
      geom.axis = ax0;
      geom.type = RABBIT_FACE;
      geom.index[0] = D - H[ax1] - 1;
      geom.index[1] = rnp[ax1] >> H[ax1];
      geom.index[2] = rnp[ax2] >> H[ax2];
      geom.index[3] = 0;
    }
    else {
      geom.axis = ax0;
      geom.type = RABBIT_EDGE;
      geom.index[0] = D - H[ax0] - 1;
      geom.index[1] = rnp[ax0] >> H[ax0];
      geom.index[2] = 0;
      geom.index[3] = 0;
    }
  }

  switch (geom.type) {
  case RABBIT_NODE:
    for (i=0; i<=1; ++i) {
      for (j=0; j<=1; ++j) {
        for (k=0; k<=1; ++k) {
          m = i*4 + j*2 + k*1;
          geom.vertices[3*m + 0] = rnp[0] + (1 << H[0]) * (i==0 ? -1 : +1);
          geom.vertices[3*m + 1] = rnp[1] + (1 << H[1]) * (j==0 ? -1 : +1);
          geom.vertices[3*m + 2] = rnp[2] + (1 << H[2]) * (k==0 ? -1 : +1);
        }
      }
    }
    break;
  case RABBIT_FACE:
    for (i=0; i<=1; ++i) {
      for (j=0; j<=1; ++j) {
	m = i*2 + j*1;
	geom.vertices[3*m + 0] = rnp[0];
	geom.vertices[3*m + 1] = rnp[1];
	geom.vertices[3*m + 2] = rnp[2];
	geom.vertices[3*m + ax1] += (1 << H[ax1]) * (i==0 ? -1 : +1);
	geom.vertices[3*m + ax2] += (1 << H[ax2]) * (j==0 ? -1 : +1);
      }
    }
    break;
  case RABBIT_EDGE:
    for (i=0; i<=1; ++i) {
      m = i*1;
      geom.vertices[3*m + 0] = rnp[0];
      geom.vertices[3*m + 1] = rnp[1];
      geom.vertices[3*m + 2] = rnp[2];
      geom.vertices[3*m + ax0] += (1 << H[ax0]) * (i==0 ? -1 : +1);
    }
    break;
  }

  return geom;
}

int rabbit_mesh_count(rabbit_mesh *M, int flags)
{
  rabbit_node *node, *tmp;
  int count = 0;

  if (flags & RABBIT_ANY) {
    return HASH_COUNT(M->nodes);
  }
  else if (flags & RABBIT_FACE) {
    return HASH_COUNT(M->faces);
  }
  else if (flags & RABBIT_EDGE) {
    return HASH_COUNT(M->edges);
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

int rabbit_mesh_merge(rabbit_mesh *M, rabbit_mesh *N)
/*
 * Merge mesh N into mesh M, preferring M on collisions. The merge only proceeds
 * if both meshes have identical config structs. Nodes, faces, and edges from N
 * are *moved* from N to M unless they already exist in M, in which case the
 * object at M is kept, and N retains ownership of it. N still needs to be freed
 * by the caller.
 */
{
  rabbit_node *node, *tmp_node, *rpl_node;
  rabbit_edge *edge, *tmp_edge, *rpl_edge;

  if (memcmp(&M->config, &N->config, sizeof(rabbit_cfg)) != 0) {
    ERR("%s", "cannot merge meshes with different config structs");
    return RABBIT_FAIL;
  }

  HASH_ITER(hh, N->nodes, node, tmp_node) {
    HASH_FIND(hh, M->nodes, node->index, 4 * sizeof(int), rpl_node);
    if (rpl_node == NULL) {
      HASH_DEL(N->nodes, node);
      HASH_ADD(hh, M->nodes, index, 4 * sizeof(int), node);
    }
  }

  HASH_ITER(hh, N->edges, edge, tmp_edge) {
    HASH_FIND(hh, M->edges, edge->vertices, 6 * sizeof(int), rpl_edge);
    if (rpl_edge == NULL) {
      HASH_DEL(N->edges, edge);
      HASH_ADD(hh, M->edges, vertices, 6 * sizeof(int), edge);
    }
  }
  return RABBIT_SUCCESS;
}

void rabbit_mesh_build(rabbit_mesh *M)
{
  rabbit_node *node, *tmp_node;
  rabbit_face *face, *tmp_face, *last_face, *existing_face;
  rabbit_edge *edge, *tmp_edge, *last_edge, *existing_edge;

  int face_rnp[3];
  int LR, a, h;
  int iter, removed;

  HASH_ITER(hh, M->nodes, node, tmp_node) {

    h = M->config.max_depth - node->index[0] - 1;

    for (LR=0; LR<=1; ++LR) {
      for (a=0; a<3; ++a) {

        face_rnp[0] = (2 * node->index[1] + 1) << h;
        face_rnp[1] = (2 * node->index[2] + 1) << h;
        face_rnp[2] = (2 * node->index[3] + 1) << h;
        face_rnp[a] = (2 * node->index[a+1] + (LR == 0 ? 0 : 2)) << h;

        HASH_FIND(hh, M->faces, face_rnp, 3 * sizeof(int), existing_face);

        if (existing_face == NULL) {

          face = (rabbit_face*) malloc(sizeof(rabbit_face));
          face->mesh = M;
          face->data = calloc(M->config.doubles_per_face, sizeof(double));
          memcpy(face->rnp, face_rnp, 3 * sizeof(int));

          HASH_ADD(hh, M->faces, rnp, 3 * sizeof(int), face);
        }
      }
    }
  }

  last_face = NULL;
  face = NULL;
  iter = 0;
  removed = 0;

  HASH_SRT(hh, M->faces, face_contiguous_compare);

  HASH_ITER(hh, M->faces, face, tmp_face) {

    int vertices[12];
    int index[3];
    int axis, depth;
    int height;
    int ax0, ax1, ax2;

    rabbit_face_geom(face, vertices, &axis, &depth);

    ax0 = axis;
    ax1 = (ax0 + 1) % 3;
    ax2 = (ax0 + 2) % 3;
    height = M->config.max_depth - depth - 1;

    index[0] = depth;
    index[1] = ((face->rnp[ax1] >> height) - 1) >> 1;
    index[2] = ((face->rnp[ax2] >> height) - 1) >> 1;

    MSG(1, "checking face %d", iter++);

    if (last_face != NULL) {
      if (face_contains(last_face, face)) {

        HASH_DEL(M->faces, last_face);
        MSG(1, "removing duplicate face %d", removed++);

        free(last_face->data);
        free(last_face);
      }
    }
    last_face = face;

  }


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
    int h = M->config.max_depth - d + 1; // height !!! fix me (+1 -> -1) !!!

    /* set the index (i,j,k) of each of this node's 8 vertices */

    vertices[0][0] = (node->index[1] + 0) << h;
    vertices[0][1] = (node->index[2] + 0) << h;
    vertices[0][2] = (node->index[3] + 0) << h;
    vertices[1][0] = (node->index[1] + 0) << h;
    vertices[1][1] = (node->index[2] + 0) << h;
    vertices[1][2] = (node->index[3] + 1) << h;
    vertices[2][0] = (node->index[1] + 0) << h;
    vertices[2][1] = (node->index[2] + 1) << h;
    vertices[2][2] = (node->index[3] + 0) << h;
    vertices[3][0] = (node->index[1] + 0) << h;
    vertices[3][1] = (node->index[2] + 1) << h;
    vertices[3][2] = (node->index[3] + 1) << h;
    vertices[4][0] = (node->index[1] + 1) << h;
    vertices[4][1] = (node->index[2] + 0) << h;
    vertices[4][2] = (node->index[3] + 0) << h;
    vertices[5][0] = (node->index[1] + 1) << h;
    vertices[5][1] = (node->index[2] + 0) << h;
    vertices[5][2] = (node->index[3] + 1) << h;
    vertices[6][0] = (node->index[1] + 1) << h;
    vertices[6][1] = (node->index[2] + 1) << h;
    vertices[6][2] = (node->index[3] + 0) << h;
    vertices[7][0] = (node->index[1] + 1) << h;
    vertices[7][1] = (node->index[2] + 1) << h;
    vertices[7][2] = (node->index[3] + 1) << h;


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
          edge->mesh = M;
          edge->data = (double*) calloc(M->config.doubles_per_edge,
                                        sizeof(double));

          memcpy(edge->vertices, edge_vertices, 6 * sizeof(int));

          HASH_ADD(hh, M->edges, vertices, 6 * sizeof(int), edge);

          MSG(2, "adding edge [%d %d %d] -> [%d %d %d] (v%d -> v%d)",
              edge->vertices[0], edge->vertices[1], edge->vertices[2],
              edge->vertices[3], edge->vertices[4], edge->vertices[5],
              v0, v1);
        }
      }
    }
  }


  last_edge = NULL;
  edge = NULL;
  iter = 0;
  removed = 0;

  HASH_SRT(hh, M->edges, edge_contiguous_compare);

  HASH_ITER(hh, M->edges, edge, tmp_edge) {

    MSG(2, "checking edge %d", iter++);

    if (last_edge != NULL) {
      if (edge_contains(last_edge, edge)) {

        HASH_DEL(M->edges, last_edge);
        MSG(2, "removing duplicate edge %d", removed++);

        free(last_edge->data);
        free(last_edge);
      }
    }
    last_edge = edge;
  }
}

void rabbit_mesh_dump(rabbit_mesh *M, char *fname)
{
  int n;
  int I[4], V[6];
  rabbit_cfg config_val;
  double node_data_val;
  double edge_data_val;
  rabbit_node *node, *tmp_node;
  rabbit_edge *edge, *tmp_edge;
  tpl_node *tn = tpl_map("S(iiii)A(i#A(f))A(i#A(f))",
                         &config_val,     // 0
                         I, 4,            // 1
                         &node_data_val,  // 2
                         V, 6,            // 3
                         &edge_data_val); // 4

  config_val = M->config;
  tpl_pack(tn, 0);

  HASH_ITER(hh, M->nodes, node, tmp_node) {
    memcpy(I, node->index, 4 * sizeof(int));

    tpl_pack(tn, 1);

    for (n=0; n<M->config.doubles_per_node; ++n) {
      node_data_val = node->data[n];
      tpl_pack(tn, 2);
    }
    tpl_pack(tn, 1); // pack data array
  }

  HASH_ITER(hh, M->edges, edge, tmp_edge) {
    memcpy(V, edge->vertices, 6 * sizeof(int));

    tpl_pack(tn, 3);

    for (n=0; n<M->config.doubles_per_edge; ++n) {
      edge_data_val = edge->data[n];
      tpl_pack(tn, 4);
    }
    tpl_pack(tn, 3); // pack data array
  }

  tpl_dump(tn, TPL_FILE, fname);
  tpl_free(tn);
}

rabbit_mesh *rabbit_mesh_load(char *fname)
{
  rabbit_mesh *M;

  int n;
  int I[4], V[6];
  rabbit_cfg config_val;
  double node_data_val;
  double edge_data_val;
  rabbit_node *node;
  rabbit_edge *edge;
  tpl_node *tn = tpl_map("S(iiii)A(i#A(f))A(i#A(f))",
                         &config_val,     // 0
                         I, 4,            // 1
                         &node_data_val,  // 2
                         V, 6,            // 3
                         &edge_data_val); // 4

  tpl_load(tn, TPL_FILE, fname);
  tpl_unpack(tn, 0);

  M = rabbit_mesh_new(config_val);

  while (tpl_unpack(tn, 1) > 0) {
    node = rabbit_mesh_putnode(M, I, RABBIT_ACTIVE);

    tpl_unpack(tn, 1); // unpack data array

    for (n=0; n<M->config.doubles_per_node; ++n) {
      tpl_unpack(tn, 2);
      node->data[n] = node_data_val;
    }
  }

  while (tpl_unpack(tn, 3) > 0) {
    edge = (rabbit_edge*) malloc(sizeof(rabbit_edge));
    edge->mesh = M;
    edge->data = (double*) calloc(M->config.doubles_per_edge,
                                  sizeof(double));
    memcpy(edge->vertices, V, 6 * sizeof(int));

    tpl_unpack(tn, 3); // unpack data array

    for (n=0; n<M->config.doubles_per_edge; ++n) {
      tpl_unpack(tn, 4);
      edge->data[n] = edge_data_val;
    }
    HASH_ADD(hh, M->edges, vertices, 6 * sizeof(int), edge);
  }

  tpl_free(tn);
  return M;
}

void rabbit_face_geom(rabbit_face *F, int vertices[12], int *axis, int *depth)
/*
 * Calculate the rational number position a face's vertices in counter-clockwise
 * order, the orientation and depth of the face
 *
 * Faces are keyed by their rational number position (rnp), three integers
 * labeling the coordinates of the face's center. The depth and orientation of
 * the face are inferred from the least significant active bit in the rnp.
 *
 */
{
  int D = F->mesh->config.max_depth;
  int H[3] = { 0, 0, 0 }; // height

  while (((F->rnp[0] >> H[0]) & 1) == 0 && H[0] < D) ++H[0];
  while (((F->rnp[1] >> H[1]) & 1) == 0 && H[1] < D) ++H[1];
  while (((F->rnp[2] >> H[2]) & 1) == 0 && H[2] < D) ++H[2];

  if (H[1] == H[2]) {

    /* x-directed face */
    if (axis) *axis = 0;
    if (depth) *depth = D - H[1] - 1;
    if (vertices) {
      vertices[ 0] = F->rnp[0];
      vertices[ 1] = F->rnp[1] - (1 << H[1]);
      vertices[ 2] = F->rnp[2] - (1 << H[1]);
      vertices[ 3] = F->rnp[0];
      vertices[ 4] = F->rnp[1] + (1 << H[1]);
      vertices[ 5] = F->rnp[2] - (1 << H[1]);
      vertices[ 6] = F->rnp[0];
      vertices[ 7] = F->rnp[1] + (1 << H[1]);
      vertices[ 8] = F->rnp[2] + (1 << H[1]);
      vertices[ 9] = F->rnp[0];
      vertices[10] = F->rnp[1] - (1 << H[1]);
      vertices[11] = F->rnp[2] + (1 << H[1]);
    }
  }
  if (H[2] == H[0]) {

    /* y-directed face */
    if (axis) *axis = 1;
    if (depth) *depth = D - H[2] - 1;
    if (vertices) {
      vertices[ 0] = F->rnp[0] - (1 << H[2]);
      vertices[ 1] = F->rnp[1];
      vertices[ 2] = F->rnp[2] - (1 << H[2]);
      vertices[ 3] = F->rnp[0] - (1 << H[2]);
      vertices[ 4] = F->rnp[1];
      vertices[ 5] = F->rnp[2] + (1 << H[2]);
      vertices[ 6] = F->rnp[0] + (1 << H[2]);
      vertices[ 7] = F->rnp[1];
      vertices[ 8] = F->rnp[2] + (1 << H[2]);
      vertices[ 9] = F->rnp[0] + (1 << H[2]);
      vertices[10] = F->rnp[1];
      vertices[11] = F->rnp[2] - (1 << H[2]);
    }
  }
  if (H[0] == H[1]) { /* z-directed face */

    /* z-directed face */
    if (axis) *axis = 2;
    if (depth) *depth = D - H[0] - 1;
    if (vertices) {
      vertices[ 0] = F->rnp[0] - (1 << H[0]);
      vertices[ 1] = F->rnp[1] - (1 << H[0]);
      vertices[ 2] = F->rnp[2];
      vertices[ 3] = F->rnp[0] + (1 << H[0]);
      vertices[ 4] = F->rnp[1] - (1 << H[0]);
      vertices[ 5] = F->rnp[2];
      vertices[ 6] = F->rnp[0] + (1 << H[0]);
      vertices[ 7] = F->rnp[1] + (1 << H[0]);
      vertices[ 8] = F->rnp[2];
      vertices[ 9] = F->rnp[0] - (1 << H[0]);
      vertices[10] = F->rnp[1] + (1 << H[0]);
      vertices[11] = F->rnp[2];
    }
  }
}

void rabbit_edge_vertices(rabbit_edge *E, int vertices[6])
{
  memcpy(vertices, E->vertices, 6 * sizeof(int));
}

int face_contiguous_compare(rabbit_face *A, rabbit_face *B)
{
  int A_index[3];
  int B_index[3];
  int A_depth;
  int B_depth;
  int A_axis;
  int B_axis;
  int A_label;
  int B_label;
  int A_height;
  int B_height;
  int ax0, ax1, ax2;

  rabbit_face_geom(A, NULL, &A_axis, &A_depth);
  rabbit_face_geom(B, NULL, &B_axis, &B_depth);

  A_height = A->mesh->config.max_depth - A_depth - 1;
  B_height = B->mesh->config.max_depth - B_depth - 1;

  /* orientation (x, y, z) - directed */
  if (A_axis != B_axis) {
    return A_axis - B_axis;
  }

  ax0 = A_axis;
  ax1 = (ax0 + 1) % 3;
  ax2 = (ax0 + 2) % 3;

  /* coordinate along face normal */
  if (A->rnp[ax0] != B->rnp[ax0]) {
    return A->rnp[ax0] - B->rnp[ax0];
  }

  A_index[0] = A_depth;
  A_index[1] = ((A->rnp[ax1] >> A_height) - 1) >> 1;
  A_index[2] = ((A->rnp[ax2] >> A_height) - 1) >> 1;

  B_index[0] = B_depth;
  B_index[1] = ((B->rnp[ax1] >> B_height) - 1) >> 1;
  B_index[2] = ((B->rnp[ax2] >> B_height) - 1) >> 1;

  /* 2d preorder label in the plane of both faces */
  A_label = preorder_label(A_index, A->mesh->config.max_depth, 2);
  B_label = preorder_label(B_index, B->mesh->config.max_depth, 2);

  return A_label - B_label;
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
    return axis_a - axis_b;
  }

  /* coordinate of next axis */
  if (A->vertices[ax1] != B->vertices[ax1]) {
    return A->vertices[ax1] - B->vertices[ax1];
  }

  /* coordinate of next axis */
  if (A->vertices[ax2] != B->vertices[ax2]) {
    return A->vertices[ax2] - B->vertices[ax2];
  }

  /* left endpoint of segment along its own axis */
  if (A->vertices[ax0] != B->vertices[ax0]) {
    return A->vertices[ax0] - B->vertices[ax0];
  }

  /* depth of segment */
  if (depth_a != depth_b) {
    return depth_a - depth_b;
  }

  return 0;
}

int face_contains(rabbit_face *A, rabbit_face *B)
/*
 * Return true if face A contains face B
 */
{
  int A_vertices[12];
  int B_vertices[12];
  int A_axis;
  int B_axis;
  int A_depth;
  int B_depth;
  int ax0, ax1, ax2;

  rabbit_face_geom(A, A_vertices, &A_axis, &A_depth);
  rabbit_face_geom(B, B_vertices, &B_axis, &B_depth);

  ax0 = A_axis;
  ax1 = (ax0 + 1) % 3;
  ax2 = (ax0 + 2) % 3;

  if (A_axis != B_axis) {
    return 0;
  }

  if (A->rnp[ax0] != B->rnp[ax0]) {
    return 0;
  }

  return (A_vertices[3*0 + ax1] <= B_vertices[3*0 + ax1] &&
          A_vertices[3*0 + ax2] <= B_vertices[3*0 + ax2] &&
          A_vertices[3*2 + ax1] >= B_vertices[3*2 + ax1] &&
          A_vertices[3*2 + ax2] >= B_vertices[3*2 + ax2]);
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

int node_preorder_compare(rabbit_node *A, rabbit_node *B)
{
  return node_preorder_label(A) - node_preorder_label(B);
}

uint64_t node_preorder_label(rabbit_node *node)
{
  return preorder_label(node->index, node->mesh->config.max_depth, 3);
}

uint64_t preorder_label(int index[4], int max_depth, int r)
/*
 * Return the order in which a given node is visited in a preorder traversal of
 * a fully fleshed out tree having arbitrary max_depth and branching ratio m=2^r
 * where r=1,2,3
 */
{
  int m = 1 << r; // branching ratio
  int d, n, h, nb, sd, Md, adding, label=0;
  uint64_t morton;

  switch (r) {
  case 1: morton = index[1]; break;
  case 2: morton = interleave_bits2(index[1], index[2]); break;
  case 3: morton = interleave_bits3(index[1], index[2], index[3]); break;
  default: morton = 0;
  }

  for (d=0; d<index[0]; ++d) {
    n = index[0] - d - 1; // bit
    h = max_depth - d - 1;
    nb = 1 << (r*n);
    sd = (morton / nb) % m;
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

uint64_t interleave_bits2(uint64_t a, uint64_t b)
/*
 * Create a 64-bit integer by interleaving the bits of the 32-bit integers a and
 * b
 */
{
  int n;
  uint64_t label = 0;

  for (n=0; n<32; ++n) {
    label += ((a >> n) & 1) << (2*n + 0);
    label += ((b >> n) & 1) << (2*n + 1);
  }

  return label;
}


#ifndef RABBIT_LIB
static void sanity_tests()
{
  /* 1d trees, m=2 */
  ASSERTEQI(tree_size_atlevel(1, 0), 1);
  ASSERTEQI(tree_size_atlevel(1, 1), 3);
  ASSERTEQI(tree_size_atlevel(1, 2), 7);

  /* 2d trees, m=4 */
  ASSERTEQI(tree_size_atlevel(2, 0), 1);
  ASSERTEQI(tree_size_atlevel(2, 1), 5);
  ASSERTEQI(tree_size_atlevel(2, 2), 21);

  /* 3d trees, m=8 */
  ASSERTEQI(tree_size_atlevel(3, 0), 1);
  ASSERTEQI(tree_size_atlevel(3, 1), 9);
  ASSERTEQI(tree_size_atlevel(3, 2), 73);

  /* does a cube have 12 edges and 6 faces? */
  if (1) {
    int I[4] = { 0, 0, 0, 0 };
    rabbit_cfg config = { 10, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_node *node = rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
    node->data[0] = 10.0;
    node->data[1] = 20.0;
    node->data[2] = 30.0;
    node->data[3] = 40.0;
    rabbit_mesh_build(mesh);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_EDGE), 12);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_FACE), 6);
    rabbit_mesh_dump(mesh, "rabbit-single.mesh");
    rabbit_mesh_del(mesh);
  }
  if (1) {
    int I[4] = { 0, 0, 0, 0 };
    rabbit_mesh *mesh = rabbit_mesh_load("rabbit-single.mesh");
    rabbit_node *node = rabbit_mesh_getnode(mesh, I);
    ASSERTEQF(node->data[0], 10.0);
    ASSERTEQF(node->data[1], 20.0);
    ASSERTEQF(node->data[2], 30.0);
    ASSERTEQF(node->data[3], 40.0);
    ASSERTEQI(mesh->config.max_depth, 10);
    ASSERTEQI(mesh->config.doubles_per_node, 4);
    ASSERTEQI(mesh->config.doubles_per_edge, 4);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_ACTIVE), 1);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_EDGE), 12);
    rabbit_mesh_del(mesh);
  }
  if (1) {
    int merge_error;
    int I0[4] = { 1, 0, 0, 0 };
    int I1[4] = { 1, 1, 0, 0 };
    rabbit_cfg config0 = { 10, 4, 4, 4 };
    rabbit_cfg config1 = { 11, 4, 4, 4 };
    rabbit_mesh *mesh0 = rabbit_mesh_new(config0);
    rabbit_mesh *mesh1 = rabbit_mesh_new(config1);
    rabbit_mesh_putnode(mesh0, I0, RABBIT_ACTIVE);
    rabbit_mesh_putnode(mesh1, I1, RABBIT_ACTIVE);
    rabbit_mesh_build(mesh0);
    rabbit_mesh_build(mesh1);
    merge_error = rabbit_mesh_merge(mesh0, mesh1);
    ASSERTEQI(merge_error, RABBIT_FAIL);
    config1.max_depth = 10;
    rabbit_mesh_del(mesh1);
    mesh1 = rabbit_mesh_new(config1);
    rabbit_mesh_putnode(mesh1, I1, RABBIT_ACTIVE);
    rabbit_mesh_build(mesh1);
    merge_error = rabbit_mesh_merge(mesh0, mesh1);
    ASSERTEQI(merge_error, RABBIT_SUCCESS);
    ASSERTEQI(rabbit_mesh_count(mesh1, RABBIT_ACTIVE), 0); // 0 overlapping nodes
    ASSERTEQI(rabbit_mesh_count(mesh1, RABBIT_EDGE), 4);   // 4 overlapping edges
    rabbit_mesh_del(mesh0);
    rabbit_mesh_del(mesh1);
  }
  if (1) {
    int D=10, d=1, i=1, j=0, k=0;
    int vertices[12];
    int axis, depth;
    rabbit_cfg config = { D, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_face F;
    F.mesh = mesh;

    /* y
     * ^
     * -------------
     * |     |  |  |
     * |     |-----|
     * |     |  |  |
     * -------------
     * |     |     |
     * |     |  o  | <--- face with rnp (3/4, 1/4, 0)
     * |     |     |      left z-face of node at depth 1: (1, 0, 0)
     * z------------> x
     */

    F.rnp[0] = (2 * i + 1) << (D - d - 1); // formula for rnp of left z-face
    F.rnp[1] = (2 * j + 1) << (D - d - 1);
    F.rnp[2] = (2 * k + 0) << (D - d - 1);

    rabbit_face_geom(&F, vertices, &axis, &depth);

    ASSERTEQI(axis, 2);
    ASSERTEQI(depth, 1);
    ASSERTEQI(vertices[0], (2 * i + 0) << (D - d - 1));
    ASSERTEQI(vertices[1], (2 * j + 0) << (D - d - 1));
    ASSERTEQI(vertices[2], (2 * k + 0) << (D - d - 1));

    ASSERTEQI(vertices[3], (2 * i + 2) << (D - d - 1));
    ASSERTEQI(vertices[4], (2 * j + 0) << (D - d - 1));
    ASSERTEQI(vertices[5], (2 * k + 0) << (D - d - 1));

    rabbit_mesh_del(mesh);
  }
  if (1) {
    int D = 3;
    rabbit_cfg config = { D, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_geom geom;
    int rnp[3];

    rnp[0] = 4;
    rnp[1] = 4;
    rnp[2] = 4;
    geom = rabbit_mesh_geom(mesh, rnp);
    ASSERTEQI(geom.type, RABBIT_NODE);
    ASSERTEQI(geom.index[0], 0);
    ASSERTEQI(geom.index[1], 0);

    rnp[0] = 2; // depth 1 node
    rnp[1] = 6;
    rnp[2] = 2;
    geom = rabbit_mesh_geom(mesh, rnp);
    ASSERTEQI(geom.type, RABBIT_NODE);
    ASSERTEQI(geom.index[0], 1);
    ASSERTEQI(geom.index[1], 0);
    ASSERTEQI(geom.index[2], 1);
    ASSERTEQI(geom.index[3], 0);

    rnp[0] = 0; // depth 1 face
    rnp[1] = 2;
    rnp[2] = 2;
    geom = rabbit_mesh_geom(mesh, rnp);
    ASSERTEQI(geom.type, RABBIT_FACE);
    ASSERTEQI(geom.axis, 0);
    ASSERTEQI(geom.index[0], 1);

    rnp[0] = 4; // x-edge at depth 0
    rnp[1] = 0;
    rnp[2] = 0;
    geom = rabbit_mesh_geom(mesh, rnp);
    ASSERTEQI(geom.type, RABBIT_EDGE);
    ASSERTEQI(geom.axis, 0);
    ASSERTEQI(geom.index[0], 0);

    rabbit_mesh_del(mesh);
  }
}

void write_meshes()
{
  /* write a uniform-depth 2d mesh */
  if (0) {
    int I[4] = { 0, 0, 0, 0 };
    int i,j;
    int d = 3;
    rabbit_cfg config = { 10, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);

    for (i=0; i<(1<<d); ++i) {
      for (j=0; j<(1<<d); ++j) {
        I[0] = d;
        I[1] = i;
        I[2] = j;
        I[3] = 0;
        rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
      }
    }
    rabbit_mesh_build(mesh);
    rabbit_mesh_dump(mesh, "rabbit-2d.mesh");
    rabbit_mesh_del(mesh);
  }
  if (1) {
    rabbit_cfg config = { 10, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    int I[4] = { 0, 0, 0, 0 };
    int i,j,k;

    for (i=0; i<2; ++i) {
      for (j=0; j<2; ++j) {
        for (k=0; k<2; ++k) {
          if (i==0 && j==0 && k==0) continue;
          I[0] = 1;
          I[1] = i;
          I[2] = j;
          I[3] = k;
          rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
        }
      }
    }

    for (i=0; i<2; ++i) {
      for (j=0; j<2; ++j) {
        for (k=0; k<2; ++k) {
          I[0] = 2;
          I[1] = i;
          I[2] = j;
          I[3] = k;
          rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
        }
      }
    }

    rabbit_mesh_build(mesh);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_FACE), 66);
    rabbit_mesh_dump(mesh, "rabbit-3d.mesh");
    rabbit_mesh_del(mesh);
  }
}

int main()
{
  sanity_tests();
  write_meshes();

  rabbit_cfg config = { 10, 4, 4, 4 };
  rabbit_mesh *mesh = rabbit_mesh_new(config);
  rabbit_node *node;
  int I[4] = { 0, 0, 0, 0 };
  int i,j,k;
  int d = 3;

  for (i=0; i<(1<<d); ++i) {
    for (j=0; j<(1<<d); ++j) {
      for (k=0; k<(1<<d); ++k) {
        I[0] = d;
        I[1] = i;
        I[2] = j;
        I[3] = k;
        node = rabbit_mesh_putnode(mesh, I, RABBIT_ACTIVE);
        node->data[0] = 0.0;
        node->data[1] = 0.1;
        node->data[2] = 0.2;
        node->data[3] = 0.3;
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
       HASH_SRT(hh, mesh->faces, face_contiguous_compare)
       );

  TIME(
       HASH_SRT(hh, mesh->edges, edge_contiguous_compare)
       );

  MSG(0, "there are %d total nodes", rabbit_mesh_count(mesh, RABBIT_ACTIVE));
  MSG(0, "there are %d total edges", rabbit_mesh_count(mesh, RABBIT_EDGE));

  rabbit_mesh_dump(mesh, "rabbit.mesh");

  rabbit_mesh_del(mesh);

  return 0;
}

#endif // RABBIT_LIB

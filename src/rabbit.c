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
static uint64_t interleave_bits2(uint64_t a, uint64_t b);
static uint64_t interleave_bits3(uint64_t a, uint64_t b, uint64_t c);
static uint64_t preorder_label(int index[4], int max_depth, int r);
static int64_t  node_preorder_compare(rabbit_node *A, rabbit_node *B);
static int64_t  face_preorder_compare(rabbit_face *A, rabbit_face *B);
static int64_t  edge_preorder_compare(rabbit_edge *A, rabbit_edge *B);
static int      face_contains(rabbit_face *A, rabbit_face *B);
static int      edge_contains(rabbit_edge *A, rabbit_edge *B);

/*
 * Return the size of a tree of max depth depth n and branching ratio m=2^r
 * tree_size_atlevel = (m^(n+1) - 1) / (m - 1)
 */
#define tree_size_atlevel(r, n) ((1 << ((r)*((n)+1))) - 1) / ((1 << (r)) - 1)

#ifdef RABBIT_USE_SYSTEM_FFS
#define FFS(i) ffs(i)
#else
static int FFS(i)
{
  int h = 0;
  while (((i >> h) & 1) == 0 && h < 32) ++h;
  return h + 1;
}
#endif // RABBIT_USE_FFS


/*
 * /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 *
 *
 *                       PUBLIC FUNCTIONS IMPLEMENTATION
 *
 *
 * \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 */

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

rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int *A, int flags)
{
  rabbit_node *node;
  int rnp[3];
  int h;

  if (flags & RABBIT_RNP) {
    rnp[0] = A[0];
    rnp[1] = A[1];
    rnp[2] = A[2];
  }
  else {
    h = M->config.max_depth - A[0] - 1;
    rnp[0] = (2 * A[1] + 1) << h;
    rnp[1] = (2 * A[2] + 1) << h;
    rnp[2] = (2 * A[3] + 1) << h;
  }

  HASH_FIND(hh, M->nodes, rnp, 3 * sizeof(int), node);

  if (node != NULL) {
    return node;
  }
  else {
    node = (rabbit_node*) malloc(sizeof(rabbit_node));
    node->flags = flags & (RABBIT_ACTIVE | RABBIT_GHOST);
    node->mesh = M;
    node->data = (double*) calloc(M->config.doubles_per_node, sizeof(double));
    memcpy(node->rnp, rnp, 3 * sizeof(int));

    HASH_ADD(hh, M->nodes, rnp, 3 * sizeof(int), node);

    return node;
  }
}

rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int *A, int flags)
{
  rabbit_node *node;

  if (flags & RABBIT_RNP) {
    HASH_FIND(hh, M->nodes, A, 3 * sizeof(int), node);
  }
  else {
    int h = M->config.max_depth - A[0] - 1;
    int rnp[3];

    rnp[0] = (2 * A[1] + 1) << h;
    rnp[1] = (2 * A[2] + 1) << h;
    rnp[2] = (2 * A[3] + 1) << h;

    HASH_FIND(hh, M->nodes, rnp, 3 * sizeof(int), node);
  }

  return node;
}

rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int *A, int flags)
{
  rabbit_node *node = rabbit_mesh_getnode(M, A, flags);

  if (node != NULL) {
    HASH_DEL(M->nodes, node);
    free(node->data);
    free(node);
  }

  return NULL;
}

rabbit_node *rabbit_mesh_containing(rabbit_mesh *M, int *A, int flags)
/*
 * Return the smallest node of which the index or rational number position (rnp)
 * A is a subset, including A if it exists
 *
 * NOTES:
 *
 *  + containment is not well-defined when the rnp is not in the open interval
 *    (0, 1)^3
 *
 *  + the rnp must represent a node, otherwise NULL will be returned
 */
{
  rabbit_node *node;
  rabbit_geom geom;
  int I[4];

  if (flags & RABBIT_RNP) {
    flags ^= (RABBIT_RNP | RABBIT_INDEX);
    geom = rabbit_mesh_geom(M, A);
    memcpy(I, geom.index, 4 * sizeof(int));
  }
  else {
    memcpy(I, A, 4 * sizeof(int));
  }

  if (I[1] < 0 || I[1] >= (1 << I[0]) ||
      I[2] < 0 || I[2] >= (1 << I[0]) || 
      I[3] < 0 || I[3] >= (1 << I[0])) {
    /* can't handle containment off-bounds */
    return rabbit_mesh_getnode(M, I, flags);
  }
  else if (I[0] < 0) {
    return NULL;
  }
  else if ((node = rabbit_mesh_getnode(M, I, flags))) {
    return node;
  }
  else {
    I[0] -= 1;
    I[1] /= 2;
    I[2] /= 2;
    I[3] /= 2;
    return rabbit_mesh_containing(M, I, flags);
  }
}

rabbit_node *rabbit_mesh_contains(rabbit_mesh *M, int *A, int flags, int *size)
/*
 * Return the subtree at index or rational number position A, and the number of
 * nodes it contains. If there are no nodes at or below the target address, NULL
 * is returned with size=0.
 */
{
  int h;
  int rnp[3];

  HASH_SRT(hh, M->nodes, node_preorder_compare);

  if (flags & RABBIT_RNP) {
    memcpy(rnp, A, 3 * sizeof(int));
    h = FFS(A[0]) - 1;
  }
  else {
    h = M->config.max_depth - A[0] - 1;
    rnp[0] = (2 * A[1] + 1) << h;
    rnp[1] = (2 * A[2] + 1) << h;
    rnp[2] = (2 * A[3] + 1) << h;
  }

  rabbit_node *iter, *head;
  rabbit_geom targ_geom = rabbit_mesh_geom(M, rnp);
  rabbit_geom iter_geom;
  uint64_t targ_label = preorder_label(targ_geom.index, M->config.max_depth, 3);
  uint64_t next_label = targ_label + tree_size_atlevel(3, h + 1);
  uint64_t iter_label;

  iter = M->nodes;

  MSG(2, "tree size at level %d is %d", h, tree_size_atlevel(3, h));

  while (iter) {
    iter_geom = rabbit_mesh_geom(M, iter->rnp);
    iter_label = preorder_label(iter_geom.index, M->config.max_depth, 3);
    if (iter_label >= targ_label) {
      break;
    }
    iter = iter->hh.next;
  }

  /* "head" node is the first one in the subtree */
  head = iter;

  *size = 0;

  while (iter) {
    iter_geom = rabbit_mesh_geom(M, iter->rnp);
    iter_label = preorder_label(iter_geom.index, M->config.max_depth, 3);

    MSG(2, "index: (%d %d) PL: %"PRIu64" upper PL: %"PRIu64,
	iter_geom.index[0], iter_geom.index[1], iter_label, next_label);

    if (iter_label >= next_label) {
      break;
    }

    *size += 1;
    iter = iter->hh.next;
  }

  if (*size == 0) {
    printf("setting to NULL, size=0\n");
    head = NULL;
  }
  return head;
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
  H[0] = rnp[0] ? FFS(rnp[0]) - 1 : D;
  H[1] = rnp[1] ? FFS(rnp[1]) - 1 : D;
  H[2] = rnp[2] ? FFS(rnp[2]) - 1 : D;

  geom.index[0] = 0;
  geom.index[1] = 0;
  geom.index[2] = 0;
  geom.index[3] = 0;

  /* if the height is the same along every axis then this is a node */
  if (H[0] == H[1] && H[1] == H[2]) {
    geom.type = RABBIT_NODE;
    geom.axis = -1; /* no orientation */
    geom.index[0] = D - H[0] - 1;
    geom.index[1] = ((rnp[0] >> H[0]) - 1) >> 1;
    geom.index[2] = ((rnp[1] >> H[0]) - 1) >> 1;
    geom.index[3] = ((rnp[2] >> H[0]) - 1) >> 1;
  }

  /*
   * otherwise it's a face or edge
   *
   * If all 3 heights are different, it's an edge and the orientation and
   * depth correspond to the smallest height.
   *
   * Put another way, it's an edge if and only if there is a unique smallest
   * height value, and otherwise it's a face. The reasoning is that getting to
   * a face from a node only involves increasing the height along one axis, so
   * the height of the other two stays the same.
   */

  else if (H[0] < H[1] && H[0] < H[2]) {
    geom.type = RABBIT_EDGE;
    geom.axis = 0;
  }
  else if (H[1] < H[2] && H[1] < H[0]) {
    geom.type = RABBIT_EDGE;
    geom.axis = 1;
  }
  else if (H[2] < H[0] && H[2] < H[1]) {
    geom.type = RABBIT_EDGE;
    geom.axis = 2;
  }
  else if (H[1] == H[2]) {
    geom.axis = 0;
    geom.type = H[0] < H[1] ? RABBIT_EDGE : RABBIT_FACE;
  }
  else if (H[2] == H[0]) {
    geom.axis = 1;
    geom.type = H[1] < H[2] ? RABBIT_EDGE : RABBIT_FACE;
  }
  else if (H[0] == H[1]) {
    geom.axis = 2;
    geom.type = H[2] < H[0] ? RABBIT_EDGE : RABBIT_FACE;
  }
  else {
    /* should never be reached, just to quiet uninitialized compiler warnings */
    geom.axis = 0;
    geom.type = RABBIT_FAIL;
  }

  if (geom.type == RABBIT_FACE) {
    ax0 = geom.axis;
    ax1 = (ax0 + 1) % 3;
    ax2 = (ax0 + 2) % 3;
    geom.index[0] = D - H[ax1] - 1;
    geom.index[1] = ((rnp[ax1] >> H[ax1]) - 1) >> 1;
    geom.index[2] = ((rnp[ax2] >> H[ax2]) - 1) >> 1;
    geom.index[3] = 0;
  }
  else if (geom.type == RABBIT_EDGE) {
    ax0 = geom.axis;
    ax1 = (ax0 + 1) % 3;
    ax2 = (ax0 + 2) % 3;
    geom.index[0] = D - H[ax0] - 1;
    geom.index[1] = ((rnp[ax0] >> H[ax0]) - 1) >> 1;
    geom.index[2] = 0;
    geom.index[3] = 0;
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

  if (flags & RABBIT_NODE) {
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
  rabbit_face *face, *tmp_face, *rpl_face;
  rabbit_edge *edge, *tmp_edge, *rpl_edge;

  if (memcmp(&M->config, &N->config, sizeof(rabbit_cfg)) != 0) {
    ERR("%s", "cannot merge meshes with different config structs");
    return RABBIT_FAIL;
  }

  HASH_ITER(hh, N->nodes, node, tmp_node) {
    HASH_FIND(hh, M->nodes, node->rnp, 3 * sizeof(int), rpl_node);
    if (rpl_node == NULL) {
      HASH_DEL(N->nodes, node);
      HASH_ADD(hh, M->nodes, rnp, 3 * sizeof(int), node);
    }
  }

  HASH_ITER(hh, N->faces, face, tmp_face) {
    HASH_FIND(hh, M->faces, face->rnp, 3 * sizeof(int), rpl_face);
    if (rpl_face == NULL) {
      HASH_DEL(N->faces, face);
      HASH_ADD(hh, M->faces, rnp, 3 * sizeof(int), face);
    }
  }

  HASH_ITER(hh, N->edges, edge, tmp_edge) {
    HASH_FIND(hh, M->edges, edge->rnp, 3 * sizeof(int), rpl_edge);
    if (rpl_edge == NULL) {
      HASH_DEL(N->edges, edge);
      HASH_ADD(hh, M->edges, rnp, 3 * sizeof(int), edge);
    }
  }
  return RABBIT_SUCCESS;
}

void rabbit_mesh_build(rabbit_mesh *M)
{
  rabbit_node *node, *tmp_node;
  rabbit_face *face, *tmp_face, *last_face;
  rabbit_edge *edge, *tmp_edge, *last_edge;

  int face_rnp[3];
  int edge_rnp[3];
  int LR, a, h;
  int iter, removed;

  HASH_ITER(hh, M->nodes, node, tmp_node) {

    rabbit_geom geom = rabbit_mesh_geom(M, node->rnp);
    h = M->config.max_depth - geom.index[0] - 1;

    for (LR=0; LR<=1; ++LR) {
      for (a=0; a<3; ++a) {

        memcpy(face_rnp, node->rnp, 3 * sizeof(int));
        face_rnp[a] += (LR == 0 ? -1 : +1) << h;

        HASH_FIND(hh, M->faces, face_rnp, 3 * sizeof(int), face);

        if (face == NULL) {

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

  HASH_SRT(hh, M->faces, face_preorder_compare);

  HASH_ITER(hh, M->faces, face, tmp_face) {

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

  HASH_ITER(hh, M->nodes, node, tmp_node) {

    int ax0, ax1, ax2;
    int i, j;

    rabbit_geom geom = rabbit_mesh_geom(M, node->rnp);
    h = M->config.max_depth - geom.index[0] - 1;

    for (a=0; a<3; ++a) {

      ax0 = (a + 0) % 3;
      ax1 = (a + 1) % 3;
      ax2 = (a + 2) % 3;

      for (i=0; i<=1; ++i) {
        for (j=0; j<=1; ++j) {

          memcpy(edge_rnp, node->rnp, 3 * sizeof(int));

          edge_rnp[ax1] += (1 << h) * (i==0 ? -1 : +1);
          edge_rnp[ax2] += (1 << h) * (j==0 ? -1 : +1);

          HASH_FIND(hh, M->edges, edge_rnp, 3 * sizeof(int), edge);

          if (edge == NULL) {

            edge = (rabbit_edge*) malloc(sizeof(rabbit_edge));
            edge->mesh = M;
            edge->data = (double*) calloc(M->config.doubles_per_edge,
                                          sizeof(double));

            memcpy(edge->rnp, edge_rnp, 3 * sizeof(int));
            HASH_ADD(hh, M->edges, rnp, 3 * sizeof(int), edge);
          }
        }
      }
    }
  }

  last_edge = NULL;
  edge = NULL;
  iter = 0;
  removed = 0;

  HASH_SRT(hh, M->edges, edge_preorder_compare);

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
  int rnp[3];
  rabbit_cfg config_val;
  double data_val;
  rabbit_node *node, *tmp_node;
  rabbit_face *face, *tmp_face;
  rabbit_edge *edge, *tmp_edge;
  tpl_node *tn = tpl_map("S(iiii)A(i#A(f))A(i#A(f))A(i#A(f))",
                         &config_val,     // 0
                         rnp, 3,          // 1 nodes
                         &data_val,       // 2
                         rnp, 3,          // 3 faces
                         &data_val,       // 4
                         rnp, 3,          // 5 edges
                         &data_val);      // 6

  config_val = M->config;
  tpl_pack(tn, 0);

  HASH_ITER(hh, M->nodes, node, tmp_node) {
    memcpy(rnp, node->rnp, 3 * sizeof(int));
    tpl_pack(tn, 1);
    for (n=0; n<M->config.doubles_per_node; ++n) {
      data_val = node->data[n];
      tpl_pack(tn, 2);
    }
    tpl_pack(tn, 1); // pack data array
  }

  HASH_ITER(hh, M->faces, face, tmp_face) {
    memcpy(rnp, face->rnp, 3 * sizeof(int));
    tpl_pack(tn, 3);
    for (n=0; n<M->config.doubles_per_face; ++n) {
      data_val = face->data[n];
      tpl_pack(tn, 4);
    }
    tpl_pack(tn, 3); // pack data array
  }

  HASH_ITER(hh, M->edges, edge, tmp_edge) {
    memcpy(rnp, edge->rnp, 3 * sizeof(int));
    tpl_pack(tn, 5);
    for (n=0; n<M->config.doubles_per_edge; ++n) {
      data_val = edge->data[n];
      tpl_pack(tn, 6);
    }
    tpl_pack(tn, 5); // pack data array
  }

  tpl_dump(tn, TPL_FILE, fname);
  tpl_free(tn);
}

rabbit_mesh *rabbit_mesh_load(char *fname)
{
  rabbit_mesh *M;

  int n;
  int rnp[3];
  rabbit_cfg config_val;
  double data_val;
  rabbit_node *node;
  rabbit_face *face;
  rabbit_edge *edge;
  tpl_node *tn = tpl_map("S(iiii)A(i#A(f))A(i#A(f))A(i#A(f))",
                         &config_val,     // 0
                         rnp, 3,          // 1 nodes
                         &data_val,       // 2
                         rnp, 3,          // 3 faces
                         &data_val,       // 4
                         rnp, 3,          // 5 edges
                         &data_val);      // 6

  tpl_load(tn, TPL_FILE, fname);
  tpl_unpack(tn, 0);

  M = rabbit_mesh_new(config_val);

  while (tpl_unpack(tn, 1) > 0) {
    node = (rabbit_node*) malloc(sizeof(rabbit_node));
    node->mesh = M;
    node->data = (double*) calloc(M->config.doubles_per_node,
                                  sizeof(double));
    node->flags = RABBIT_ACTIVE;
    memcpy(node->rnp, rnp, 3 * sizeof(int));
    tpl_unpack(tn, 1); // unpack data array
    for (n=0; n<M->config.doubles_per_node; ++n) {
      tpl_unpack(tn, 2);
      node->data[n] = data_val;
    }
    HASH_ADD(hh, M->nodes, rnp, 3 * sizeof(int), node);
  }

  while (tpl_unpack(tn, 3) > 0) {
    face = (rabbit_face*) malloc(sizeof(rabbit_face));
    face->mesh = M;
    face->data = (double*) calloc(M->config.doubles_per_face,
                                  sizeof(double));
    memcpy(face->rnp, rnp, 3 * sizeof(int));
    tpl_unpack(tn, 3); // unpack data array
    for (n=0; n<M->config.doubles_per_face; ++n) {
      tpl_unpack(tn, 4);
      face->data[n] = data_val;
    }
    HASH_ADD(hh, M->faces, rnp, 3 * sizeof(int), face);
  }

  while (tpl_unpack(tn, 5) > 0) {
    edge = (rabbit_edge*) malloc(sizeof(rabbit_edge));
    edge->mesh = M;
    edge->data = (double*) calloc(M->config.doubles_per_edge,
                                  sizeof(double));
    memcpy(edge->rnp, rnp, 3 * sizeof(int));
    tpl_unpack(tn, 5); // unpack data array
    for (n=0; n<M->config.doubles_per_edge; ++n) {
      tpl_unpack(tn, 6);
      edge->data[n] = data_val;
    }
    HASH_ADD(hh, M->edges, rnp, 3 * sizeof(int), edge);
  }

  tpl_free(tn);
  return M;
}

/*
 * /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 *
 *
 *                           INTERNAL FUNCTIONS
 *
 *
 * \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 */

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

uint64_t preorder_label(int index[4], int max_depth, int r)
/*
 * Return the order in which a given node is visited in a preorder traversal of
 * a fully fleshed out tree having arbitrary max_depth and branching ratio m=2^r
 * where r=1,2,3
 */
{
  int m = 1 << r; // branching ratio
  int d, n, h;
  uint64_t nb, sd, Md, adding, morton, label=0;

  switch (r) {
  case 1: morton = index[1]; break;
  case 2: morton = interleave_bits2(index[1], index[2]); break;
  case 3: morton = interleave_bits3(index[1], index[2], index[3]); break;
  default: morton = 0; /* to quiet compiler warnings */
  }

  for (d=0; d<index[0]; ++d) {
    n = index[0] - d - 1; // height above
    h = max_depth - d - 1; // total height
    nb = 1 << (r*n);
    sd = (morton / nb) % m;
    Md = tree_size_atlevel(r, h);
    adding = sd * Md + 1;
    label += adding;
  }

  return label;
}

int64_t node_preorder_compare(rabbit_node *A, rabbit_node *B)
{
  rabbit_geom A_geom = rabbit_mesh_geom(A->mesh, A->rnp);
  rabbit_geom B_geom = rabbit_mesh_geom(B->mesh, B->rnp);
  return (preorder_label(A_geom.index, A->mesh->config.max_depth, 3) -
          preorder_label(B_geom.index, B->mesh->config.max_depth, 3));
}

int64_t face_preorder_compare(rabbit_face *A, rabbit_face *B)
{
  uint64_t A_label;
  uint64_t B_label;
  int ax0, ax1, ax2;

  rabbit_geom A_geom = rabbit_mesh_geom(A->mesh, A->rnp);
  rabbit_geom B_geom = rabbit_mesh_geom(B->mesh, B->rnp);

  /* orientation (x, y, z) - directed */
  if (A_geom.axis != B_geom.axis) {
    return A_geom.axis - B_geom.axis;
  }

  ax0 = A_geom.axis;
  ax1 = (ax0 + 1) % 3;
  ax2 = (ax0 + 2) % 3;

  /* coordinate along face normal */
  if (A->rnp[ax0] != B->rnp[ax0]) {
    return A->rnp[ax0] - B->rnp[ax0];
  }

  /* 2d preorder label in the plane of both faces */
  A_label = preorder_label(A_geom.index, A->mesh->config.max_depth, 2);
  B_label = preorder_label(B_geom.index, B->mesh->config.max_depth, 2);

  return A_label - B_label;
}

int64_t edge_preorder_compare(rabbit_edge *A, rabbit_edge *B)
{
  uint64_t A_label;
  uint64_t B_label;
  int ax0, ax1, ax2;

  rabbit_geom A_geom = rabbit_mesh_geom(A->mesh, A->rnp);
  rabbit_geom B_geom = rabbit_mesh_geom(B->mesh, B->rnp);

  /* orientation (x, y, z) - directed */
  if (A_geom.axis != B_geom.axis) {
    return A_geom.axis - B_geom.axis;
  }

  ax0 = A_geom.axis;
  ax1 = (ax0 + 1) % 3;
  ax2 = (ax0 + 2) % 3;

  /* coordinate of next axis */
  if (A->rnp[ax1] != B->rnp[ax1]) {
    return A->rnp[ax1] - B->rnp[ax1];
  }

  /* coordinate of next axis */
  if (A->rnp[ax2] != B->rnp[ax2]) {
    return A->rnp[ax2] - B->rnp[ax2];
  }

  /* 1d preorder label on the axis of both edges */
  A_label = preorder_label(A_geom.index, A->mesh->config.max_depth, 1);
  B_label = preorder_label(B_geom.index, B->mesh->config.max_depth, 1);

  return A_label - B_label;
}

int face_contains(rabbit_face *A, rabbit_face *B)
/*
 * Return true if face A contains face B
 */
{
  int ax0, ax1, ax2;
  rabbit_geom A_geom = rabbit_mesh_geom(A->mesh, A->rnp);
  rabbit_geom B_geom = rabbit_mesh_geom(B->mesh, B->rnp);

  ax0 = A_geom.axis;
  ax1 = (ax0 + 1) % 3;
  ax2 = (ax0 + 2) % 3;

  if (A_geom.axis != B_geom.axis) {
    return 0;
  }

  if (A->rnp[ax0] != B->rnp[ax0]) {
    return 0;
  }

  return (A_geom.vertices[3*0 + ax1] <= B_geom.vertices[3*0 + ax1] &&
          A_geom.vertices[3*0 + ax2] <= B_geom.vertices[3*0 + ax2] &&
          A_geom.vertices[3*3 + ax1] >= B_geom.vertices[3*3 + ax1] &&
          A_geom.vertices[3*3 + ax2] >= B_geom.vertices[3*3 + ax2]);
}

int edge_contains(rabbit_edge *A, rabbit_edge *B)
/*
 * Return true if edge A contains edge B
 */
{
  int ax0, ax1, ax2;
  rabbit_geom A_geom = rabbit_mesh_geom(A->mesh, A->rnp);
  rabbit_geom B_geom = rabbit_mesh_geom(B->mesh, B->rnp);

  ax0 = A_geom.axis;
  ax1 = (ax0 + 1) % 3; // only used if axis_a == axis_b
  ax2 = (ax0 + 2) % 3;

  /* different orientation? */
  if (A_geom.axis != B_geom.axis) {
    return 0;
  }

  /* not co-linear? */
  if (A->rnp[ax1] != B->rnp[ax1]) {
    return 0;
  }
  if (A->rnp[ax2] != B->rnp[ax2]) {
    return 0;
  }

  return (A_geom.vertices[ax0+0] <= B_geom.vertices[ax0+0] &&
          A_geom.vertices[ax0+3] >= B_geom.vertices[ax0+3]);
}


/*
 * /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 *
 *
 *                                  TESTS
 *
 *
 * \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 */
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
    int R[3] = { 0, 0, 0 };
    rabbit_mesh *mesh = rabbit_mesh_load("rabbit-single.mesh");
    rabbit_node *node = rabbit_mesh_getnode(mesh, I, RABBIT_INDEX);
    rabbit_geom geom;
    ASSERTEQF(node->data[0], 10.0);
    ASSERTEQF(node->data[1], 20.0);
    ASSERTEQF(node->data[2], 30.0);
    ASSERTEQF(node->data[3], 40.0);
    ASSERTEQI(mesh->config.max_depth, 10);
    ASSERTEQI(mesh->config.doubles_per_node, 4);
    ASSERTEQI(mesh->config.doubles_per_edge, 4);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_ACTIVE), 1);
    ASSERTEQI(rabbit_mesh_count(mesh, RABBIT_EDGE), 12);

    I[0] = 1;
    node = rabbit_mesh_containing(mesh, I, RABBIT_INDEX);
    geom = rabbit_mesh_geom(mesh, node->rnp);
    ASSERTEQI(geom.index[0], 0);

    R[0] = 511;
    R[1] = 511;
    R[2] = 511;

    node = rabbit_mesh_containing(mesh, R, RABBIT_RNP);
    geom = rabbit_mesh_geom(mesh, node->rnp);
    ASSERTEQI(geom.index[0], 0);

    I[0] = 0;
    rabbit_mesh_delnode(mesh, I, RABBIT_INDEX);
    node = rabbit_mesh_containing(mesh, R, RABBIT_RNP);
    ASSERTEQI((node == NULL), 1);

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
    int *vertices;
    int axis, depth;
    int rnp[3];
    rabbit_cfg config = { D, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_geom geom;

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

    rnp[0] = (2 * i + 1) << (D - d - 1); // formula for rnp of left z-face
    rnp[1] = (2 * j + 1) << (D - d - 1);
    rnp[2] = (2 * k + 0) << (D - d - 1);

    geom = rabbit_mesh_geom(mesh, rnp);
    axis = geom.axis;
    depth = geom.index[0];
    vertices = geom.vertices;

    /* Demonstrates ordering of vertices around face. For CCW ordering, take
       vertices 0, 2, 3, 1 */
    ASSERTEQI(axis, 2);
    ASSERTEQI(depth, 1);
    ASSERTEQI(vertices[ 0], (2 * i + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 1], (2 * j + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 2], (2 * k + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 3], (2 * i + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 4], (2 * j + 2) << (D - d - 1));
    ASSERTEQI(vertices[ 5], (2 * k + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 6], (2 * i + 2) << (D - d - 1));
    ASSERTEQI(vertices[ 7], (2 * j + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 8], (2 * k + 0) << (D - d - 1));
    ASSERTEQI(vertices[ 9], (2 * i + 2) << (D - d - 1));
    ASSERTEQI(vertices[10], (2 * j + 2) << (D - d - 1));
    ASSERTEQI(vertices[11], (2 * k + 0) << (D - d - 1));

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
  if (1) {
    /* ensure that preorder labels respect periodic mesh */
    int D = 4;
    rabbit_cfg config = { D, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_node *n1, *n2;
    rabbit_geom g1, g2;
    int index[4];

    index[0] = 2;
    index[1] = -1;
    index[2] = -1;
    index[3] = -1;
    n1 = rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    g1 = rabbit_mesh_geom(mesh, n1->rnp);

    index[0] = 2;
    index[1] = 3;
    index[2] = 3;
    index[3] = 3;
    n2 = rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    g2 = rabbit_mesh_geom(mesh, n2->rnp);

    int L1 = preorder_label(g1.index, D, 3);
    int L2 = preorder_label(g2.index, D, 3);
    ASSERTEQI(L1, L2);

    rabbit_mesh_build(mesh);
    rabbit_mesh_dump(mesh, "rabbit-guard.mesh");
    rabbit_mesh_del(mesh);
  }
  if (1) {
    /*
     * This test uses the following tree:
     * -----------------------------------------
     *           x
     *      x         x
     *   x     o   o     x
     *  o o             o o
     * -----------------------------------------
     */
    int D = 4;
    rabbit_cfg config = { D, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_node *node;
    int index[4] = { 0, 0, 0, 0 };
    int size;

    index[0] = 3; index[1] = 0;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    index[0] = 3; index[1] = 1;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    index[0] = 2; index[1] = 1;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    index[0] = 2; index[1] = 2;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    index[0] = 3; index[1] = 6;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);
    index[0] = 3; index[1] = 7;
    rabbit_mesh_putnode(mesh, index, RABBIT_ACTIVE);


    /* check all the nodes having subtree size 1 */

    index[0] = 3; index[1] = 0;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);
    index[0] = 3; index[1] = 1;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);
    index[0] = 2; index[1] = 1;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);
    index[0] = 2; index[1] = 2;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);
    index[0] = 3; index[1] = 6;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);
    index[0] = 3; index[1] = 7;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 1);

    /* check nodes with subtree size 2 */
    index[0] = 2; index[1] = 0;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 2);
    index[0] = 2; index[1] = 3;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 2);

    /* check nodes with subtree size 3 */
    index[0] = 1; index[1] = 0;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 3);
    index[0] = 1; index[1] = 1;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 3);
    
    /* check root node, with subtree size 6 */
    index[0] = 0; index[1] = 0;
    node = rabbit_mesh_contains(mesh, index, RABBIT_INDEX, &size);
    ASSERTEQI((node != NULL), 1);
    ASSERTEQI(size, 6);
    
    rabbit_mesh_del(mesh);
  }
}

void write_meshes()
{
  /* write a uniform-depth 2d mesh */
  if (1) {
    int I[4] = { 0, 0, 0, 0 };
    int i,j;
    int d = 6;
    rabbit_cfg config = { 12, 8, 0, 0 };
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
    TIME( rabbit_mesh_build(mesh) );
    TIME( rabbit_mesh_dump(mesh, "rabbit-2d.mesh") );
    TIME( rabbit_mesh_del(mesh) );
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
  if (1) {
    rabbit_cfg config = { 10, 4, 4, 4 };
    rabbit_mesh *mesh = rabbit_mesh_new(config);
    rabbit_node *node;
    int I[4] = { 0, 0, 0, 0 };
    int i,j,k;
    int d = 5;

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

    TIME( rabbit_mesh_build(mesh) );
    TIME( HASH_SRT(hh, mesh->nodes, node_preorder_compare) );
    TIME( HASH_SRT(hh, mesh->faces, face_preorder_compare) );
    TIME( HASH_SRT(hh, mesh->edges, edge_preorder_compare) );

    MSG(0, "there are %d total nodes", rabbit_mesh_count(mesh, RABBIT_ACTIVE));
    MSG(0, "there are %d total faces", rabbit_mesh_count(mesh, RABBIT_FACE));
    MSG(0, "there are %d total edges", rabbit_mesh_count(mesh, RABBIT_EDGE));

    rabbit_mesh_dump(mesh, "rabbit-deep.mesh");
    rabbit_mesh_del(mesh);
  }
}

int main()
{
  sanity_tests();
  write_meshes();
  return 0;
}

#endif // RABBIT_LIB

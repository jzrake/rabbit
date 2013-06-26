#include <stdio.h>
#include <stdlib.h>
#include "zmesh.h"

static void realloc_faces(struct zmesh *M, int max_faces)
{
  M->faces = (struct zface*) realloc(M->faces, max_faces*sizeof(struct zface));
  M->max_faces = max_faces;
}

static void push_face(struct zmesh *M, struct ztree *nL, struct ztree *nR)
{
  struct zface F;
  if (M->num_faces >= M->max_faces) {
    realloc_faces(M, M->max_faces * 2);
  }
  F.cellL = nL;
  F.cellR = nR;
  M->faces[M->num_faces] = F;
  M->num_faces += 1;
}

struct zmesh *zmesh_new(struct ztree *tree)
/*
 * Create a new mesh from the given tree
 */
{
  struct zmesh *M = (struct zmesh*) malloc(sizeof(struct zmesh));
  M->tree = tree;
  M->faces = NULL;
  M->num_faces = 0;
  realloc_faces(M, 1);
  return M;
}

void zmesh_del(struct zmesh *M)
{
  free(M->faces);
  free(M);
}

void zmesh_build_faces(struct zmesh *M)
/*
 * This function builds faces, 2-dimensional areas that are shared by exactly
 * two volumes. A given node in the tree is responsible for building its right
 * x, y, and z-directed faces when
 *
 * (1) It is a leaf
 * (2) The volume to its right is at the same refinement level or finer
 *
 * That node is responsible for building its left x, y, and z-directed faces
 * when
 *
 * (1) It is a leaf
 * (2) The volume to its left is at a strictly coarser refinement level
 *
 * These requirements ensure that faces which divide volumes of different
 * refinement levels are built for the finer volume.
 */
{
  zmesh_clear_faces(M);
  switch (ztree_rank(M->tree)) {
  case 1: {
    struct ztree *it = NULL, *nL, *nR;
    while ((it = ztree_next_leaf(M->tree, it))) {
      nL = ztree_travel1(it, 0, -1);
      nR = ztree_travel1(it, 0, +1);
      if (nL) {
	/* the volume to the left has a strictly coarser refinement level */
	if (ztree_depth(nL) < ztree_depth(it) && ztree_isleaf(nL)) {
	  push_face(M, nL, it);
	}
      }
      if (nR) {
	/* the volume to the right has a coarser or equal refinement level */
	if (ztree_depth(nR) <= ztree_depth(it) && ztree_isleaf(nR)) {
	  push_face(M, it, nR);
	}
      }
    }
    break;
  }
  case 2: /* NOT YET IMPLEMENTED */ break;
  case 3: /* NOT YET IMPLEMENTED */ break;
  }
}

void zmesh_clear_faces(struct zmesh *M)
{
  M->num_faces = 0;
  realloc_faces(M, 1);
}

int zmesh_num_faces(const struct zmesh *M)
{
  return M->num_faces;
}

struct zface zmesh_get_face(const struct zmesh *M, int n)
{
  return M->faces[n];
}

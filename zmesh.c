#include <stdio.h>
#include <stdlib.h>
#include "zmesh.h"

struct zmesh *zmesh_new(struct ztree *tree)
/*
 * Create a new mesh from the given tree
 */
{
  struct zmesh *M = (struct zmesh *) malloc(sizeof(struct zmesh));
  M->tree = tree;
  M->faces = NULL;
  M->num_faces = 0;
  return M;
}

void zmesh_build_faces(struct zmesh *M)
{

}

void zmesh_del(struct zmesh *M)
{
  free(M->faces);
  free(M);
}

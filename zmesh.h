#ifndef _ZMESH_INCLUDE_
#define _ZMESH_INCLUDE_
#include "ztree.h"

struct zmesh
{
  struct ztree *tree;
  struct zface *faces;
  int num_faces;
} ;
struct zface
{
  struct ztree *cellL;
  struct ztree *cellR;
  struct zface *next;
} ;

struct zmesh *zmesh_new(struct ztree *tree);
void zmesh_del(struct zmesh *M);
void zmesh_build_faces(struct zmesh *M);

#endif // _ZMESH_INCLUDE_

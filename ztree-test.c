#include <stdio.h>
#include <assert.h>
#include "ztree.h"
#include "zmesh.h"

#define nohup 0
#define asserteq(E, v)printf("%s == %d : %d\n",#E,v,E); assert(E==v||nohup);

static void output_tree(const struct ztree *T, const char *fname);
static void output_faces(const struct zmesh *M, const char *fname);

int test_tree()
{
  struct ztree *tree = ztree_new(3, sizeof(double));
  struct ztree *it = NULL;

  asserteq(ztree_rank(tree), 3);
  asserteq(ztree_descendant_node_count(tree), 0);
  asserteq(ztree_descendant_leaf_count(tree), 0);

  ztree_split(tree);
  asserteq(ztree_descendant_node_count(tree), 8);
  asserteq(ztree_descendant_leaf_count(tree), 8);

  ztree_prune(tree);
  asserteq(ztree_descendant_node_count(tree), 0);
  asserteq(ztree_descendant_leaf_count(tree), 0);

  ztree_split(tree);
  ztree_split(tree);
  asserteq(ztree_descendant_node_count(tree), 72);
  asserteq(ztree_descendant_leaf_count(tree), 64);

  int I[3] = { 0, 1, 2 };
  struct ztree *child = ztree_travel(tree, 2, I);

  asserteq(ztree_depth(child), 2);
  asserteq(ztree_index(child, 0), 0);
  asserteq(ztree_index(child, 1), 1);
  asserteq(ztree_index(child, 2), 2);

  ztree_prune(tree);
  ztree_split(tree);
  ztree_split(tree);
  int node_count = 0; // will include root node (=72 + 1)
  it = NULL;
  while ((it = ztree_next(tree, it))) ++node_count;
  asserteq(ztree_descendant_node_count(tree), node_count - 1);
  ztree_del(tree);

  /* Test for more complex 1D traversal */
  tree = ztree_new(1, 0);
  ztree_prune(tree);
  ztree_splitn(tree, 5);
  ztree_split(ztree_travel1(tree, 3, 2));
  ztree_split(ztree_travel1(tree, 3, 5));
  it = ztree_travel1(tree, 6, 16);
  asserteq(ztree_index(it, 0), 16);
  asserteq(ztree_depth(it), 6);
  ztree_travel1(it, 0, -1);
  ztree_del(tree);

  return 0;
}

int test_mesh()
{
  struct ztree *tree = ztree_new(1, sizeof(double));
  struct zmesh *mesh = zmesh_new(tree);

  ztree_splitn(tree, 5);
  zmesh_build_faces(mesh);
  asserteq(zmesh_num_faces(mesh), 31);

  ztree_prune(tree);
  zmesh_build_faces(mesh);
  asserteq(zmesh_num_faces(mesh), 0);

  ztree_prune(tree);
  ztree_splitn(tree, 5);
  ztree_split(ztree_travel1(tree, 3, 2));
  ztree_split(ztree_travel1(tree, 3, 5));
  zmesh_build_faces(mesh);
  asserteq(zmesh_num_faces(mesh), 39);

  output_tree(tree, "tree.dat");
  output_faces(mesh, "faces.dat");

  zmesh_del(mesh);
  ztree_del(tree);

  return 0;
}

void output_tree(const struct ztree *T, const char *fname)
{
  FILE *outf = fopen(fname, "w");
  struct ztree *it = NULL;
  while ((it = ztree_next(T, it))) {
    double x = (ztree_index(it, 0) + 0.5) / (1 << ztree_depth(it));
    fprintf(outf, "%f %d\n", x, ztree_depth(it));
  }
}

void output_faces(const struct zmesh *M, const char *fname)
{
  FILE *outf = fopen(fname, "w");
  int n;
  for (n=0; n<zmesh_num_faces(M); ++n) {
    struct zface F = zmesh_get_face(M, n);
    double xL = (ztree_index(F.cellL, 0) + 0.5) / (1 << ztree_depth(F.cellL));
    double xR = (ztree_index(F.cellR, 0) + 0.5) / (1 << ztree_depth(F.cellR));
    double z = 0.5 * (ztree_depth(F.cellL) + ztree_depth(F.cellR));
    fprintf(outf, "%f %f\n", 0.5 * (xL + xR), z);
  }
}

int main()
{
  test_tree();
  test_mesh();
  return 0;
}

#include <stdio.h>
#include <assert.h>
#include "ztree.h"
#include "zmesh.h"

#define nohup 0
#define asserteq(E, v)printf("%s == %d : %d\n",#E,v,E); assert(E==v||nohup);

//static void output_tree(const struct ztree *T, const char *fname);
//static void output_faces(const struct zmesh *M, const char *fname);

int test_tree()
{
  struct ztree *tree = ztree_new(1, 0);
  struct ztree *it;
  int n, I[1];
  int nodes = 0;

  ztree_branch(tree);
  asserteq(ztree_count(tree, ZTREE_ROOT), 1);
  asserteq(ztree_count(tree, ZTREE_STUB), 2);
  asserteq(ztree_count(tree, ZTREE_LEAF), 0);

  ztree_splitn(tree, 2);

  asserteq(ztree_count(tree, ZTREE_ROOT), 1);
  asserteq(ztree_count(tree, ZTREE_STUB), 0);
  asserteq(ztree_count(tree, ZTREE_LEAF), 4);

  it = ztree_travel1(tree, 2, 0);
  asserteq(ztree_depth(it), 2);
  asserteq(ztree_index(it, 0), 0);

  I[0] = 0;
  it = ztree_require_node(it, 2, I);
  asserteq(ztree_depth(it), 4);
  asserteq(ztree_status(ztree_require_node(tree, -1, I)), ZTREE_STUB);

  it = NULL;
  while ((it = ztree_next(tree, it))) {
    nodes += 1;
  }
  asserteq(ztree_count(tree, ZTREE_NODE), nodes);

  ztree_prune(tree);
  asserteq(ztree_count(tree, ZTREE_NODE), 1);
  for (n=0; n < 1<<12; ++n) {
    it = ztree_require_node(tree, 12, &n);
  }
  asserteq(ztree_count(tree, ZTREE_LEAF), 1<<12);

  ztree_del(tree);
  return 0;
}

int test_mesh()
{
  return 0;
}

int main(int argc, char **argv)
{
  test_tree();
  test_mesh();
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
  fclose(outf);
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
  fclose(outf);
}

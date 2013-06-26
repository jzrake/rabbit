#include <stdio.h>
#include <assert.h>
#include "ztree.h"
#include "zmesh.h"

#define nohup 1
#define asserteq(E, v)printf("%s == %s : %d\n",#E,#v,E); assert(E==v||nohup);

int test_tree()
{
  struct ztree *tree = ztree_new(3, sizeof(double));

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
  struct ztree *it = NULL;
  int node_count = 0;
  while ((it = ztree_next(tree, it))) {
    ++node_count;
  }

  asserteq(ztree_descendant_node_count(tree), node_count);
  ztree_del(tree);

  return 0;
}

int test_mesh()
{
  struct ztree *tree = ztree_new(1, sizeof(double));
  struct zmesh *mesh = zmesh_new(tree);
  ztree_split(tree);
  ztree_split(tree);
  ztree_split(tree);
  ztree_split(tree);
  zmesh_build_faces(mesh);
  zmesh_del(mesh);
  ztree_del(tree);

  return 0;
}

int main()
{
  test_tree();
  test_mesh();
  return 0;
}

#include <stdlib.h>
#define _ZTREE_PRIVATE_
#include "ztree.h"

#define IS_ROOT (T->parent == NULL)
#define IS_LEAF (T->children == NULL)
#define MAX_RANK 32

struct ztree *ztree_new(unsigned int rank)
/*
 * Create a new tree node with given rank
 */
{
  struct ztree *T = (struct ztree *) malloc(sizeof(struct ztree));
  T->id = 0;
  T->bytes = 0;
  T->rank = rank;
  T->data = NULL;
  T->parent = NULL;
  T->children = NULL;
  return T;
}

void ztree_del(struct ztree *T)
/*
 * Free T and all descendant nodes recursively
 */
{
  ztree_prune(T);
  free(T->data);
  free(T);
}

void ztree_prune(struct ztree *T)
/*
 * Remove descendant nodes recursively, return silently if T is a leaf node
 */
{
  unsigned int n;
  if (T->children == NULL) return;
  for (n=0; n < 1<<T->rank; ++n) {
    ztree_del(T->children[n]);
  }
  free(T->children);
  T->children = NULL;
}

void ztree_split(struct ztree *T)
/*
 * Populate child nodes immediately below T
 */
{
  unsigned int n;
  if (IS_LEAF) {
    T->children = (struct ztree **) malloc((1<<T->rank) * sizeof(struct ztree *));
    for (n=0; n < 1<<T->rank; ++n) {
      T->children[n] = ztree_new(T->rank);
      T->children[n]->parent = T;
      T->children[n]->id = n;
    }
  }
  else {
    for (n=0; n < 1<<T->rank; ++n) {
      ztree_split(T->children[n]);
    }
  }
}

int ztree_index(const struct ztree *T, int axis)
/*
 * Return the node index relative to the root node along the given axis
 */
{
  if (IS_ROOT || axis >= T->rank) {
    return 0;
  }
  else {
    return ztree_index(T->parent, axis) * 2 + (1 & (T->id >> axis));
  }
}

int ztree_id(const struct ztree *T)
{
  return T->id;
}

int ztree_descendant_node_count(const struct ztree *T)
/*
 * Return the total number of descendants below the present node
 */
{
  unsigned int i;
  unsigned int n = 1 << T->rank;
  if (IS_LEAF) return 0;
  for (i=0; i < 1<<T->rank; ++i) {
    n += ztree_descendant_node_count(T->children[i]);
  }
  return n;
}

int ztree_descendant_leaf_count(const struct ztree *T)
/*
 * Return the total number of leaf nodes below the present node
 */
{
  unsigned int i;
  unsigned int n = 0;
  if (IS_ROOT && IS_LEAF) return 0;
  if (IS_LEAF) return 1;
  for (i=0; i < 1<<T->rank; ++i) {
    n += ztree_descendant_leaf_count(T->children[i]);
  }
  return n;
}

int ztree_rank(const struct ztree *T)
{
  return T->rank;
}

int ztree_depth(const struct ztree *T)
{
  return IS_ROOT ? 0 : ztree_depth(T->parent) + 1;
}

struct ztree *ztree_parent(const struct ztree *T)
{
  return T->parent;
}

struct ztree *ztree_next(const struct ztree *T, const struct ztree *P)
/*
 * Traverse the tree depth first. Leaf nodes return their next sibling, except
 * the last one which returns its uncle.
 */
{
  /* first time through */
  if (P == NULL) {
    return (struct ztree*) T;
  }
  /* do nothing if the starting node is a leaf. */
  else if (IS_LEAF) {
    return NULL;
  }
  /* visit your children if you have any, starting with the youngest */
  else if (P->children != NULL) {
    return P->children[0];
  }
  else {
    /* move on to the next sibling if you are not the oldest child */
    if (P->id != (1<<P->rank)-1) {
      return P->parent->children[P->id + 1];
    }
    /* if you are the oldest sibling go to your closest ancestor who is not */
    else {
      while (P->id == (1<<P->rank)-1) {
        P = P->parent;
	/* if any ancestor along the way is the starting node, stop there */
        if (P == T) {
	  return NULL;
        }
      }
      /* this child is not the oldest; go on to its next sibling */
      return P->parent->children[P->id + 1];
    }
  }
}

struct ztree *ztree_travel(const struct ztree *T, int depth, const int *I0)
/*
 * Travel 'depth' levels farther from the root, and across I0 nodes at the
 * target depth. If any component of I0 is negative, or larger than 2^depth,
 * then the target node is not a descendant of this node. In that case, climb to
 * the closest common ancestor and descend the tree to the target node. Depth
 * and breadth traversals are taken relative to the present node, not the root
 * node.
 */
{
  struct ztree *p = (struct ztree *) T;
  int d, pid;
  int target_out_of_range = 0;
  int I[MAX_RANK];
  for (d=0; d<T->rank; ++d) I[d] = I0[d];

  if (depth < 0) {
    if (IS_ROOT) {
      return NULL;
    }
    else {
      return ztree_travel(T->parent, depth + 1, I0);
    }
  }

  for (d=0; d<T->rank; ++d) {
    /*
     * If an index is negative or greater than 2^depth, the target node is not
     * our descendant.
     */
    target_out_of_range += (0 > I[d] || I[d] >= (1 << depth));
  }
  if (target_out_of_range) {
    if (IS_ROOT) {
      return NULL;
    }
    else {
      /*
       * The correction below adds 2^depth to the index offset when the present
       * node has id=1 along the d-axis. This accounts for the fact that
       * deferring to the parent either walks backward in index or stays the
       * same, but never forward.
       */
      for (d=0; d<T->rank; ++d) {
        I[d] += ((1 & (T->id >> d)) << depth);
      }
      return ztree_travel(T->parent, depth + 1, I);
    }
  }
  /*
   * If we made it this far, then the target node is a descendant of the present
   * node. Locating it is easy; just step left if id=0 along the d-axis, and
   * right if id=1, until we reach the desired depth.
   */
  while (depth--) {
    if (IS_LEAF) {
      return p; /* silently return the closest ancestor */
    }
    else {
      pid = 0;
      for (d=0; d<T->rank; ++d) {
        pid += ((1 & (I[d] >> depth)) << d);
      }
      p = p->children[pid];
    }
  }
  return p;
}

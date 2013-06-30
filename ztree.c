/*
 * FILE: ztree.c
 *
 * AUTHOR: Jonathan Zrake
 *
 * DESCRIPTION:
 *
 * Tree are built using a "bottom-up" tree scheme. Leaves may be added to any
 * level below a given node by specifying the depth and index. Intermediate
 * nodes are created as necessary.
 *
 * Traversals may return branches or leafs. If the target node is not present,
 * but the search terminates at a leaf, then that leaf, the "closest ancestor"
 * is returned. If the search terminates at a stub, then NULL is returned.
 *
 * When a node is "branched" its child pointer is allocated but none of the
 * children are created. It is then no longer a leaf but a branch. Child nodes
 * are only created through the add_leaf function.
 *
 * Deleting a node first deletes all of its descendants, then clears itself from
 * its parent's child list. If all children have been deleted from the parent's
 * child list then the parent's child list is cleared and the parent returns to
 * being a leaf.
 *
 * A stub is a where a branch contains a child node set to NULL. The stub
 * implies a tree boundary, which may indicate that the tree continues on a
 * remote processor. If the stub data is required, a remote query can be issued
 * to which all processors will respond with a code indicating one of the
 * following:
 *
 * 0: I don't have that node at all
 * 1: I don't have that node, but it resides below one of my leaves
 * 2: I have that node, and it's a branch
 * 3: I have that node, and it's a leaf
 *
 * Leaf status is inferred: (T->children == NULL) <=> IS_LEAF
 *
 * Root status is inferred: (T->parent == NULL) <=> IS_ROOT
 *
 * Stub status is inferred: (T->children[id] == NULL)
 *
 * Dead branch status is inferred: (T->children[n] == NULL for all n)
 */

#include <stdlib.h>
#include <string.h>
#define _ZTREE_PRIVATE_
#include "ztree.h"

#define IS_ROOT (T->parent == NULL)
#define IS_LEAF (T->children == NULL)

static int create_intermediate_nodes = 0;

struct ztree *ztree_new(unsigned int rank, unsigned int bytes)
/*
 * Create a new tree node with given rank
 */
{
  struct ztree *T = (struct ztree*) malloc(sizeof(struct ztree));
  T->id = 0;
  T->bytes = bytes;
  T->data = malloc(bytes);
  T->rank = rank;
  T->parent = NULL;
  T->children = NULL;
  return T;
}

void ztree_del(struct ztree *T)
/*
 * Free T and all descendant nodes recursively, clear self from parent's list of
 * children
 */
{
  if (T == NULL) {
    return;
  }
  if (!IS_ROOT) {
    T->parent->children[T->id] = NULL;
  }
  ztree_prune(T);
  free(T->data);
  free(T);
}

enum ztree_node_status ztree_status(const struct ztree *T)
/*
 * Return the status of a given node. Every node has exactly one of the possible
 * status flags: stub, leaf, dead, or branch.
 */
{
  if (T == NULL) {
    return ZTREE_STUB; // it's a stub
  }
  else if (T->children == NULL) {
    return ZTREE_LEAF;
  }
  else {
    int n, dead=1;
    for (n=0; n < 1<<T->rank; ++n) {
      dead *= (T->children[n] == NULL);
    }
    if (dead) {
      return ZTREE_DEAD;
    }
    else {
      return ZTREE_BRANCH;
    }
  }
}

void ztree_prune(struct ztree *T)
/*
 * Remove descendant nodes recursively, return silently if T is a leaf node
 */
{
  unsigned int n;
  if (IS_LEAF) return;
  for (n=0; n < 1<<T->rank; ++n) {
    ztree_del(T->children[n]);
  }
  free(T->children);
  T->children = NULL; // T is now a leaf
}

void ztree_split(struct ztree *T)
/*
 * Populate all child nodes immediately below T, turning child stubs
 * into leaves. Recurse into each pre-existing child node.
 */
{
  unsigned int n;
  if (IS_LEAF) {
    ztree_branch(T);
  }
  for (n=0; n < 1<<T->rank; ++n) {
    if (T->children[n] == NULL) {
      T->children[n] = ztree_new(T->rank, T->bytes);
      T->children[n]->parent = T;
      T->children[n]->id = n;
    }
    else {
      ztree_split(T->children[n]);
    }
  }
}

void ztree_splitn(struct ztree *T, int n)
/*
 * Split the tree n times (convenience function)
 */
{
  while (n--) ztree_split(T);
}

void ztree_branch(struct ztree *T)
/*
 * If T is a leaf, then allocate its list of children, leaving them all as
 * stubs. If T is a branch then recurse to and branch all descendant leaf nodes.
 */
{
  unsigned int n;
  if (IS_LEAF) {
    T->children = (struct ztree **) malloc((1<<T->rank) * sizeof(struct ztree *));
    for (n=0; n < 1<<T->rank; ++n) {
      T->children[n] = NULL;
    }
  }
  else {
    for (n=0; n < 1<<T->rank; ++n) {
      if (T->children[n] != NULL) {
        ztree_branch(T->children[n]);
      }
    }
  }
}

void ztree_get_data_buffer(const struct ztree *T, void **buffer)
{
  *buffer = T->data;
}

void ztree_address(const struct ztree *T, struct zaddress *A)
/*
 * Return the vector index I and the depth relative to the root node
 */
{
  int d;
  for (d=0; d<T->rank; ++d) {
    A->index[d] = ztree_index(T, d);
  }
  for (d=T->rank; d<ZTREE_MAX_RANK; ++d) {
    A->index[d] = 0;
  }
  A->depth = ztree_depth(T);
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

int ztree_count(const struct ztree *T, enum ztree_node_status type)
{
  unsigned int i;
  unsigned int n = 0;
  switch (type) {
  case ZTREE_BRANCH: n += !IS_LEAF; break;
  case ZTREE_LEAF: n += IS_LEAF; break;
  case ZTREE_NODE: n += 1; break;
  case ZTREE_ROOT: n += IS_ROOT; break;
  case ZTREE_DEAD: n += ztree_status(T) == ZTREE_DEAD; break;
  case ZTREE_STUB:
    if (!IS_LEAF) {
      for (i=0; i < 1<<T->rank; ++i) {
        n += (T->children[i] == NULL);
      }
    }
    break;
  }
  if (!IS_LEAF) {
    for (i=0; i < 1<<T->rank; ++i) {
      if (T->children[i] != NULL) {
        n += ztree_count(T->children[i], type);
      }
    }
  }
  return n;
}

int ztree_depth(const struct ztree *T)
{
  return IS_ROOT ? 0 : ztree_depth(T->parent) + 1;
}

int ztree_rank(const struct ztree *T)
{
  return T->rank;
}

int ztree_isleaf(const struct ztree *T)
{
  return IS_LEAF;
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
  int n;
  /* first time through */
  if (P == NULL) {
    return (struct ztree*) T;
  }
  /* do nothing if T is a leaf */
  else if (IS_LEAF) {
    return NULL;
  }
  /* visit P's first child if it's a branch */
  else if (P->children != NULL) {
    for (n=0; n < 1<<T->rank; ++n) {
      if (P->children[n]) {
        return P->children[n];
      }
    }
  }
  /* P is a leaf or a dead branch; go on to the next sibling if there is one */
  while (ztree_next_sibling(P) == NULL) {
    P = P->parent;
    /* if any ancestor along the way is the starting node, stop there */
    if (P == T) {
      return NULL;
    }
  }
  return ztree_next_sibling(P);
}

struct ztree *ztree_next_leaf(const struct ztree *T, const struct ztree *P)
{
  while (1) {
    P = ztree_next(T, P);
    if (P == NULL) return (struct ztree*) P;
    if (ztree_isleaf(P)) return (struct ztree*) P;
  }
}

struct ztree *ztree_require_node(struct ztree *T, int depth, const int *I)
/*
 * Create (or locate) and return the target node relative to T, creating
 * intermediate branches as necessary. If the target node is out-of-bounds, or
 * already exists then NULL is returned.
 */
{
  struct ztree *leaf;
  create_intermediate_nodes = 1;
  leaf = ztree_travel(T, depth, I);
  create_intermediate_nodes = 0;
  return leaf;
}

struct ztree *ztree_next_sibling(const struct ztree *T)
/*
 * Return the next child of T's parent node, and NULL if T is the last one
 */
{
  int n;
  for (n=T->id + 1; n < 1<<T->rank; ++n) {
    if (T->parent->children[n]) {
      return T->parent->children[n];
    }
  }
  return NULL;
}

struct ztree *ztree_travel(const struct ztree *T, int depth, const int *I0)
/*
 * Travel 'depth' levels farther from the root, and across I0 nodes at the
 * target depth. If any component of I0 is negative, or larger than 2^depth,
 * then the target node is not a descendant of this node. In that case, climb to
 * the closest common ancestor and descend the tree to the target node. Depth
 * and breadth traversals are taken relative to the present node, not the root
 * node. If the target node is off the tree then return NULL. If the target node
 * is on the tree but not occupied then its nearest ancestor is returned.
 */
{
  struct ztree *p = (struct ztree *) T;
  int d, pid;
  int target_out_of_range = 0;
  int I[ZTREE_MAX_RANK];
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
    if (p->children == NULL) {
      if (create_intermediate_nodes) {
        ztree_branch(p);
      }
      else {
        /* return the closest ancestor of the target node if it does not exist */
        return p;
      }
    }
    pid = 0;
    for (d=0; d<T->rank; ++d) {
      pid += ((1 & (I[d] >> depth)) << d);
    }
    if (p->children[pid] == NULL) {
      if (create_intermediate_nodes) {
        p->children[pid] = ztree_new(T->rank, T->bytes);
        p->children[pid]->parent = p;
        p->children[pid]->id = pid;
      }
      else {
        return NULL;
      }
    }
    p = p->children[pid];
  }
  return p;
}

struct ztree *ztree_travel1(const struct ztree *T, int d, int i)
{
  int I[1] = { i };
  return ztree_travel(T, d, I);
}

struct ztree *ztree_travel2(const struct ztree *T, int d, int i, int j)
{
  int I[2] = { i, j };
  return ztree_travel(T, d, I);
}

struct ztree *ztree_travel3(const struct ztree *T, int d, int i, int j, int k)
{
  int I[3] = { i, j, k };
  return ztree_travel(T, d, I);
}

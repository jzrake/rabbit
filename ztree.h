#ifndef _ZTREE_INCLUDE_
#define _ZTREE_INCLUDE_

struct ztree;

struct ztree *ztree_new(unsigned int rank);
void ztree_del(struct ztree *T);
void ztree_prune(struct ztree *T);
void ztree_split(struct ztree *T);
int ztree_index(const struct ztree *T, int axis);
int ztree_id(const struct ztree *T);
int ztree_descendant_node_count(const struct ztree *T);
int ztree_descendant_leaf_count(const struct ztree *T);
int ztree_depth(const struct ztree *T);
int ztree_rank(const struct ztree *T);
struct ztree *ztree_parent(const struct ztree *T);
struct ztree *ztree_next(const struct ztree *T, const struct ztree *S);
struct ztree *ztree_travel(const struct ztree *T, int depth, const int *I0);

#ifdef _ZTREE_PRIVATE_
struct ztree
{
  void *data; // data buffer
  unsigned int bytes; // number of bytes for data
  unsigned int rank; // number of dimensions in tree
  unsigned int id; // child index relative to parent node
  struct ztree *parent;
  struct ztree **children;
} ;
#endif // _ZTREE_PRIVATE_
#endif // _ZTREE_INCLUDE_

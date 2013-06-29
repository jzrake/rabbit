#ifndef _ZTREE_INCLUDE_
#define _ZTREE_INCLUDE_

#define ZTREE_MAX_RANK 32

enum ztree_node_type {
  ZTREE_ROOT,
  ZTREE_BRANCH,
  ZTREE_LEAF,
  ZTREE_STUB,
  ZTREE_GHOST,
} ;

struct ztree;
struct zaddress;

struct ztree *ztree_new(unsigned int rank, unsigned int bytes);
void ztree_del(struct ztree *T);

void ztree_prune(struct ztree *T);
void ztree_split(struct ztree *T);
void ztree_splitn(struct ztree *T, int n);
void ztree_get_data_buffer(const struct ztree *T, void **buffer);
void ztree_address(const struct ztree *T, struct zaddress *A);
int ztree_index(const struct ztree *T, int axis);
int ztree_id(const struct ztree *T);
int ztree_descendant_node_count(const struct ztree *T);
int ztree_descendant_leaf_count(const struct ztree *T);
int ztree_depth(const struct ztree *T);
int ztree_rank(const struct ztree *T);
int ztree_isleaf(const struct ztree *T);

struct ztree *ztree_parent(const struct ztree *T);
struct ztree *ztree_next(const struct ztree *T, const struct ztree *S);
struct ztree *ztree_next_leaf(const struct ztree *T, const struct ztree *P);
struct ztree *ztree_add_leaf(struct ztree *T, int depth, const int *I);
struct ztree *ztree_travel(const struct ztree *T, int depth, const int *I0);
struct ztree *ztree_travel1(const struct ztree *T, int d, int i);
struct ztree *ztree_travel2(const struct ztree *T, int d, int i, int j);
struct ztree *ztree_travel3(const struct ztree *T, int d, int i, int j, int k);

struct zaddress
{
  int depth;
  int index[ZTREE_MAX_RANK];
} ;

#ifdef _ZTREE_PRIVATE_
struct ztree
{
  enum ztree_node_type type;
  void *data; // data buffer
  unsigned int bytes; // number of bytes for data
  unsigned int rank; // number of dimensions in tree
  unsigned int id; // child index relative to parent node
  struct ztree *parent;
  struct ztree **children;
} ;
#endif // _ZTREE_PRIVATE_
#endif // _ZTREE_INCLUDE_

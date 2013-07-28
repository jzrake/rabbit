#ifndef RABBIT_H
#define RABBIT_H

#include <inttypes.h> /* uint64_t */

/* use ffs call to identify least significant bits, may not be portable */
#define RABBIT_USE_SYSTEM_FFS

/* ------------------------------------------------------------------
 * RABBIT PUBLIC API
 * --------------------------------------------------------------- */
#define RABBIT_SUCCESS  0
#define RABBIT_FAIL     (1 << 0)
#define RABBIT_ANY      (1 << 1)
#define RABBIT_ACTIVE   (1 << 2)
#define RABBIT_GHOST    (1 << 3)
#define RABBIT_FORCE    (1 << 4)
#define RABBIT_NODE     (1 << 5)
#define RABBIT_FACE     (1 << 6)
#define RABBIT_EDGE     (1 << 7)
#define RABBIT_RNP      (1 << 8)
#define RABBIT_INDEX    (1 << 9)

typedef struct rabbit_mesh rabbit_mesh;
typedef struct rabbit_node rabbit_node;
typedef struct rabbit_face rabbit_face;
typedef struct rabbit_edge rabbit_edge;
typedef struct
{
  int max_depth;
  int doubles_per_node;
  int doubles_per_face;
  int doubles_per_edge;
} rabbit_cfg;
typedef struct
{
  int index[4];            /* index (depth, i, j, k) */
  int vertices[24];        /* 3 int's x 8(node), 4(face), 2(edge) */
  int type;                /* node, face, or edge */
  int axis;                /* edge direction, or face normal */
  uint64_t preorder_label; /* 3d for node, 2d for face, 1d for edge */
} rabbit_geom;

void         rabbit_mesh_del(rabbit_mesh *M);
rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg);
rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int *A, int flags);
rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int *A, int flags);
rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int *A, int flags);
rabbit_node *rabbit_mesh_containing(rabbit_mesh *M, int *A, int flags);
rabbit_geom  rabbit_mesh_geom(rabbit_mesh *M, int rnp[3]);
int          rabbit_mesh_count(rabbit_mesh *M, int flags);
int          rabbit_mesh_merge(rabbit_mesh *M, rabbit_mesh *N);
void         rabbit_mesh_build(rabbit_mesh *M);
void         rabbit_mesh_dump(rabbit_mesh *M, char *fname);
rabbit_mesh *rabbit_mesh_load(char *fname);


/* ------------------------------------------------------------------
 * RABBIT INTERNALS
 * --------------------------------------------------------------- */
#ifdef RABBIT_INTERNAL

/* debug level options */
#define ALWAYS 3
#define SOMETIMES 2
#define ALMOST_NEVER 1
#define NEVER 0

/* print messages how often? */
#define PRINT_MESSAGES ALMOST_NEVER

#include <stdio.h>
#define MSG(level, format, ...) do {            \
    if (level < PRINT_MESSAGES) {               \
      fprintf(stdout, "[%s]$ ",  __FUNCTION__); \
      fprintf(stdout, format, __VA_ARGS__);     \
      fprintf(stdout, "\n");                    \
    }                                           \
  } while (0)                                   \

#define ERR(format, ...) do {				\
    fprintf(stderr, "[ERROR:%s]$ ",  __FUNCTION__);	\
    fprintf(stderr, format, __VA_ARGS__);		\
    fprintf(stderr, "\n");				\
  } while (0)						\

#include <time.h>
#define TIME(cmd) do {                                  \
    clock_t start = clock();                            \
    cmd;                                                \
    printf("[%s]$ %s took %5.4f ms\n", __FUNCTION__,    \
           #cmd, 1e3*(clock() - start)/CLOCKS_PER_SEC); \
  } while (0)


/* assertion macro */
#include <assert.h>
#define NOHUP 0 // continue even if an assertion fails
#define ASSERT_I "[assertion:%s]$ %s == %s : %d\n"
#define ASSERT_F "[assertion:%s]$ %s == %s : %f\n"
#define ASSERTEQI(E,v)printf(ASSERT_I,__FUNCTION__,#E,#v,E);assert(E==v||NOHUP);
#define ASSERTEQF(E,v)printf(ASSERT_F,__FUNCTION__,#E,#v,E);assert(E==v||NOHUP);


/* ----------------------------------
 * Hash functions available in uthash
 * ----------------------------------
 * JEN Jenkins (default)
 * BER Bernstein
 * SAX Shift-Add-Xor
 * OAT One-at-a-time
 * FNV Fowler/Noll/Vo
 * SFH Paul Hsieh
 * MUR MurmurHash v3 (see note)
 */
#define HASH_FUNCTION HASH_JEN
#include "uthash.h"
#include "tpl.h" /* tpl file format header */


/* ----------------------------------
 * Data structure definitions
 * ----------------------------------
 */
struct rabbit_mesh {
  rabbit_cfg config;
  rabbit_node *nodes;
  rabbit_face *faces;
  rabbit_edge *edges;
} ;

struct rabbit_node {
  int rnp[3];
  int flags;
  double *data;
  rabbit_mesh *mesh;
  UT_hash_handle hh;
} ;

struct rabbit_face {
  int rnp[3]; // rational number position
  double *data;
  rabbit_mesh *mesh;
  UT_hash_handle hh;
} ;

struct rabbit_edge {
  int rnp[3]; // rational number position
  double *data;
  rabbit_mesh *mesh;
  UT_hash_handle hh;
} ;

#endif // RABBIT_INTERNAL
#endif // RABBIT_H

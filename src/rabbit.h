#ifndef RABBIT_H
#define RABBIT_H

/* ------------------------------------------------------------------
 * RABBIT PUBLIC API
 * --------------------------------------------------------------- */
#define RABBIT_ANY      (1 << 0)
#define RABBIT_ACTIVE   (1 << 1)
#define RABBIT_GHOST    (1 << 2)
#define RABBIT_FORCE    (1 << 3)
#define RABBIT_EDGE     (1 << 4)


typedef struct rabbit_mesh rabbit_mesh; // opaque
typedef struct rabbit_node rabbit_node;
typedef struct rabbit_face rabbit_face;
typedef struct rabbit_edge rabbit_edge;
typedef struct // fully defined
{
  int max_depth;
  int doubles_per_node;
  int doubles_per_edge;
} rabbit_cfg;


void         rabbit_mesh_del(rabbit_mesh *M);
rabbit_mesh *rabbit_mesh_new(rabbit_cfg cfg);
rabbit_node *rabbit_mesh_putnode(rabbit_mesh *M, int index[4], int flags);
rabbit_node *rabbit_mesh_getnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_delnode(rabbit_mesh *M, int index[4]);
rabbit_node *rabbit_mesh_containing(rabbit_mesh *M, int index[4]);
int          rabbit_mesh_count(rabbit_mesh *M, int flags);
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
      fprintf(stderr, "[%s]$ ",  __FUNCTION__); \
      fprintf(stderr, format, __VA_ARGS__);     \
      fprintf(stderr, "\n");                    \
    }                                           \
  } while (0)                                   \

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
#define ASSERT_MSG "[assertion:%s]$ %s == %d : %d\n"
#define ASSERTEQ(E,v)printf(ASSERT_MSG,__FUNCTION__,#E,v,E);assert(E==v||NOHUP);


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
#define HASH_FUNCTION HASH_SAX
#include "uthash.h"
#include "tpl.h" /* tpl file format header */


/* ----------------------------------
 * Data structure definitions
 * ----------------------------------
 */
struct rabbit_mesh {
  rabbit_cfg config;
  rabbit_node *nodes;
  rabbit_edge *edges;
} ;

struct rabbit_node {
  int index[4]; // (depth, i, j, k)
  int flags;
  double *data;
  rabbit_mesh *mesh;
  UT_hash_handle hh;
} ;

struct rabbit_face {
  rabbit_node nodes[2];
} ;

struct rabbit_edge {
  int vertices[6];
  double *data;
  UT_hash_handle hh;
} ;

#endif // RABBIT_INTERNAL
#endif // RABBIT_H

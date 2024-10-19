#include "vector.h"

/* a curve is a set of pt2s produced by a function taking a real ranging
 * between 0 and 1. curves can be linked together sequentially. */
typedef struct curve {
    pt2(*func)(struct curve*, real);
    struct curve* neighbors[2];
} curve;

/* a surface is a set of pt3s produced by a function taking a pt2 ranging
 * between <0,0> and <1,1>. surfaces can be linked in up to 4 directions,
 * fewer if some edges are of 0 length. */
typedef struct surface {
     pt3(*func)(struct surface*, pt2);
    struct surface* neighbors[2][2];
} surface;

surface** surface_neighbor(surface* sface, int dim_0_or_1, int edge_0_or_1)
{
    return &sface->neighbors[dim_0_or_1][edge_0_or_1];
}

/* a surface group collects linked surfaces together. */
typedef struct surface_group
{
    surface ** surfaces;
    unsigned int n_surfaces;
    void (*free)(struct surface_group *);
} surface_group;

/* a kind of surface defined by a position and a tangent plane. */
typedef struct oriented_surface
{
    struct surface surface;
    pt3 position;
    pt3 axes[2];
} oriented_surface;

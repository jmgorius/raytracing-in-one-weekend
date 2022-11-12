#ifndef INCLUDED_AABB_H
#define INCLUDED_AABB_H

#include "point3.h"
#include "ray.h"

#include <stdbool.h>

typedef struct AABB {
  Point3 min;
  Point3 max;
} AABB;

bool aabb_hit(const AABB *aabb, Ray r, double t_min, double t_max);
AABB aabb_surrounding_box(const AABB *first, const AABB *second);

#endif /* INCLUDED_AABB_H */

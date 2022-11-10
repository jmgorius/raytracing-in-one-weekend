#ifndef INCLUDED_POINT3_H
#define INCLUDED_POINT3_H

#include "vec3.h"

typedef struct point3 {
  double x, y, z;
} point3;

point3 point3_add(point3 p, vec3 v);
vec3 point3_sub(point3 p1, point3 p2);

#endif /* INCLUDED_POINT3_H */

#ifndef INCLUDED_RAY_H
#define INCLUDED_RAY_H

#include "point3.h"
#include "vec3.h"

typedef struct ray {
  point3 origin;
  vec3 direction;
} ray;

point3 ray_at(ray r, double t);

#endif /* INCLUDED_RAY_H */

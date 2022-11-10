#ifndef INCLUDED_RAY_H
#define INCLUDED_RAY_H

#include "point3.h"
#include "vec3.h"

typedef struct Ray {
  Point3 origin;
  Vec3 direction;
} Ray;

Point3 ray_at(Ray r, double t);

#endif /* INCLUDED_RAY_H */

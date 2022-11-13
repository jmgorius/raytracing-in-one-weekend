#ifndef INCLUDED_POINT3_H
#define INCLUDED_POINT3_H

#include "vec3.h"

typedef struct Point3 {
  double x, y, z;
} Point3;

Point3 point3_add(Point3 p, Vec3 v);
Vec3 point3_sub(Point3 p1, Point3 p2);

Point3 point3_random_in_range(double min, double max);

#endif /* INCLUDED_POINT3_H */

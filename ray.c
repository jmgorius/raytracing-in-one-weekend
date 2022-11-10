#include "ray.h"
#include "point3.h"

Point3 ray_at(Ray r, double t) {
  return point3_add(r.origin, vec3_mul(t, r.direction));
}

#include "ray.h"
#include "point3.h"

point3 ray_at(ray r, double t) {
  return point3_add(r.origin, vec3_mul(t, r.direction));
}

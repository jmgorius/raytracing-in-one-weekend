#include "aabb.h"

#include <math.h>

bool aabb_hit(const AABB *aabb, Ray r, double t_min, double t_max) {
  /* X */
  {
    double inv_d = 1.0 / r.direction.x;
    double t0 = (aabb->min.x - r.origin.x) * inv_d;
    double t1 = (aabb->max.x - r.origin.x) * inv_d;
    if (inv_d < 0.0) {
      double tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    t_min = t0 > t_min ? t0 : t_min;
    t_max = t1 < t_max ? t1 : t_max;
    if (t_max <= t_min)
      return false;
  }
  /* Y */
  {
    double inv_d = 1.0 / r.direction.y;
    double t0 = (aabb->min.y - r.origin.y) * inv_d;
    double t1 = (aabb->max.y - r.origin.y) * inv_d;
    if (inv_d < 0.0) {
      double tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    t_min = t0 > t_min ? t0 : t_min;
    t_max = t1 < t_max ? t1 : t_max;
    if (t_max <= t_min)
      return false;
  }
  /* Z */
  {
    double inv_d = 1.0 / r.direction.z;
    double t0 = (aabb->min.z - r.origin.z) * inv_d;
    double t1 = (aabb->max.z - r.origin.z) * inv_d;
    if (inv_d < 0.0) {
      double tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    t_min = t0 > t_min ? t0 : t_min;
    t_max = t1 < t_max ? t1 : t_max;
    if (t_max <= t_min)
      return false;
  }

  return true;
}

AABB aabb_surrounding_box(const AABB *first, const AABB *second) {
  return (AABB){
      .min = {fmin(first->min.x, second->min.x),
              fmin(first->min.y, second->min.y),
              fmin(first->min.z, second->min.z)},
      .max = {fmax(first->max.x, second->max.x),
              fmax(first->max.y, second->max.y),
              fmax(first->max.z, second->max.z)},
  };
}

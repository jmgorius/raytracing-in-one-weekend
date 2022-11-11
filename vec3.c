#include "vec3.h"
#include "utils.h"

#include <math.h>

Vec3 vec3_neg(Vec3 v) { return (Vec3){-v.x, -v.y, -v.z}; }

Vec3 vec3_add(Vec3 v1, Vec3 v2) {
  return (Vec3){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

Vec3 vec3_sub(Vec3 v1, Vec3 v2) {
  return (Vec3){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

Vec3 vec3_mul(double c, Vec3 v) { return (Vec3){c * v.x, c * v.y, c * v.z}; }

Vec3 vec3_div(Vec3 v, double c) { return (Vec3){v.x / c, v.y / c, v.z / c}; }

double vec3_length(Vec3 v) { return sqrt(vec3_length2(v)); }

double vec3_length2(Vec3 v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

Vec3 vec3_normalize(Vec3 v) { return vec3_div(v, vec3_length(v)); }

double vec3_dot(Vec3 v1, Vec3 v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec3 vec3_cross(Vec3 v1, Vec3 v2) {
  return (Vec3){
      v1.y * v2.z - v1.z * v2.y,
      v1.z * v2.x - v1.x * v2.z,
      v1.x * v2.y - v1.y * v2.x,
  };
}

Vec3 vec3_random(void) {
  return (Vec3){
      random_double(),
      random_double(),
      random_double(),
  };
}

Vec3 vec3_random_in_range(double min, double max) {
  return (Vec3){
      random_double_in_range(min, max),
      random_double_in_range(min, max),
      random_double_in_range(min, max),
  };
}

Vec3 vec3_random_in_unit_sphere(void) {
  while (1) {
    Vec3 result = vec3_random_in_range(-1, 1);
    if (vec3_length2(result) >= 1)
      continue;
    return result;
  }
}

Vec3 vec3_random_unit_vector(void) {
  return vec3_normalize(vec3_random_in_unit_sphere());
}

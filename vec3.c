#include "vec3.h"

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

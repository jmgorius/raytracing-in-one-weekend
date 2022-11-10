#include "vec3.h"

#include <math.h>

vec3 vec3_neg(vec3 v) { return (vec3){-v.x, -v.y, -v.z}; }

vec3 vec3_add(vec3 v1, vec3 v2) {
  return (vec3){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

vec3 vec3_sub(vec3 v1, vec3 v2) {
  return (vec3){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

vec3 vec3_mul(double c, vec3 v) { return (vec3){c * v.x, c * v.y, c * v.z}; }

vec3 vec3_div(vec3 v, double c) { return (vec3){v.x / c, v.y / c, v.z / c}; }

double vec3_length(vec3 v) { return sqrt(vec3_length2(v)); }

double vec3_length2(vec3 v) { return v.x * v.x + v.y * v.y + v.z * v.z; }

vec3 vec3_normalize(vec3 v) { return vec3_div(v, vec3_length(v)); }

double vec3_dot(vec3 v1, vec3 v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec3 vec3_cross(vec3 v1, vec3 v2) {
  return (vec3){
      v1.y * v2.z - v1.z * v2.y,
      v1.z * v2.x - v1.x * v2.z,
      v1.x * v2.y - v1.y * v2.x,
  };
}

#ifndef INCLUDED_VEC3_H
#define INCLUDED_VEC3_H

typedef struct vec3 {
  double x, y, z;
} vec3;

vec3 vec3_neg(vec3 v);
vec3 vec3_add(vec3 v1, vec3 v2);
vec3 vec3_sub(vec3 v1, vec3 v2);
vec3 vec3_mul(double c, vec3 v);
vec3 vec3_div(vec3 v, double c);

double vec3_length(vec3 v);
double vec3_length2(vec3 v);

vec3 vec3_normalize(vec3 v);

double vec3_dot(vec3 v1, vec3 v2);
vec3 vec3_cross(vec3 v1, vec3 v2);

#endif /* INCLUDED_VEC3_H */

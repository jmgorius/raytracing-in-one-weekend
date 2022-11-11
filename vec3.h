#ifndef INCLUDED_VEC3_H
#define INCLUDED_VEC3_H

typedef struct Vec3 {
  double x, y, z;
} Vec3;

Vec3 vec3_neg(Vec3 v);
Vec3 vec3_add(Vec3 v1, Vec3 v2);
Vec3 vec3_sub(Vec3 v1, Vec3 v2);
Vec3 vec3_mul(double c, Vec3 v);
Vec3 vec3_div(Vec3 v, double c);

double vec3_length(Vec3 v);
double vec3_length2(Vec3 v);

Vec3 vec3_normalize(Vec3 v);

double vec3_dot(Vec3 v1, Vec3 v2);
Vec3 vec3_cross(Vec3 v1, Vec3 v2);

Vec3 vec3_random(void);
Vec3 vec3_random_in_range(double min, double max);
Vec3 vec3_random_in_unit_sphere(void);
Vec3 vec3_random_unit_vector(void);

#endif /* INCLUDED_VEC3_H */

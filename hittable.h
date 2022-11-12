#ifndef INCLUDED_HITTABLE_H
#define INCLUDED_HITTABLE_H

#include "arena.h"
#include "material.h"
#include "point3.h"
#include "ray.h"
#include "vec3.h"

#include <stdbool.h>
#include <stddef.h>

typedef struct HitRecord {
  const Material *material;
  Point3 p;
  Vec3 normal;
  double t;
  bool front_face;
} HitRecord;

void hit_record_set_face_normal(HitRecord *record, Ray r, Vec3 outward_normal);

typedef enum HittableType {
  HITTABLE_LIST,
  HITTABLE_SPHERE,
} HittableType;

typedef struct Hittable {
  HittableType type;
} Hittable;

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record);

typedef struct HittableList {
  HittableType type;
  const Hittable **objects;
  size_t size;
  size_t capacity;
} HittableList;

void hittable_list_add(HittableList *list, const Hittable *hittable,
                       Arena *arena);
bool hittable_list_hit(const HittableList *list, Ray r, double t_min,
                       double t_max, HitRecord *record);

typedef struct Sphere {
  HittableType type;
  const Material *material;
  Point3 center;
  double radius;
} Sphere;

bool sphere_hit(const Sphere *sphere, Ray r, double t_min, double t_max,
                HitRecord *record);

#endif /* INCLUDED_HITTABLE_H */

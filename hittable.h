#ifndef INCLUDED_HITTABLE_H
#define INCLUDED_HITTABLE_H

#include "aabb.h"
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
  HITTABLE_MOVING_SPHERE,
  HITTABLE_BVH_NODE,
} HittableType;

typedef struct Hittable {
  HittableType type;
} Hittable;

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record);
bool hittable_bounding_box(const Hittable *hittable, double time_start,
                           double time_end, AABB *bounding_box);

typedef struct HittableList {
  HittableType type;
  const Hittable **objects;
  size_t size;
  size_t capacity;
} HittableList;

HittableList *hittable_list_create(Arena *arena);
void hittable_list_add(HittableList *list, const Hittable *hittable,
                       Arena *arena);
bool hittable_list_hit(const HittableList *list, Ray r, double t_min,
                       double t_max, HitRecord *record);
bool hittable_list_bounding_box(const HittableList *list, double time_start,
                                double time_end, AABB *bounding_box);

typedef struct Sphere {
  HittableType type;
  const Material *material;
  Point3 center;
  double radius;
} Sphere;

Sphere *sphere_create(Point3 center, double radius, const Material *material,
                      Arena *arena);
bool sphere_hit(const Sphere *sphere, Ray r, double t_min, double t_max,
                HitRecord *record);
bool sphere_bounding_box(const Sphere *sphere, double time_start,
                         double time_end, AABB *bounding_box);

typedef struct MovingSphere {
  HittableType type;
  const Material *material;
  Point3 center_start, center_end;
  double start, end;
  double radius;
} MovingSphere;

MovingSphere *moving_sphere_create(Point3 center_start, Point3 center_end,
                                   double start, double end, double radius,
                                   const Material *material, Arena *arena);
Point3 moving_sphere_center(const MovingSphere *sphere, double t);
bool moving_sphere_hit(const MovingSphere *sphere, Ray r, double t_min,
                       double t_max, HitRecord *record);
bool moving_sphere_bounding_box(const MovingSphere *sphere, double time_start,
                                double time_end, AABB *bounding_box);

typedef struct BVHNode {
  HittableType type;
  const Hittable *left;
  const Hittable *right;
  AABB box;
} BVHNode;

BVHNode *bvh_node_create(const Hittable **objects, size_t start, size_t end,
                         double time_start, double time_end, Arena *arena);
bool bvh_node_hit(const BVHNode *node, Ray r, double t_min, double t_max,
                  HitRecord *record);
bool bvh_node_bounding_box(const BVHNode *node, double time_start,
                           double time_end, AABB *bounding_box);

#endif /* INCLUDED_HITTABLE_H */

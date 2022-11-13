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
  double u, v;
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

typedef struct Hittable Hittable;

typedef struct HittableList {
  const Hittable **objects;
  size_t size;
  size_t capacity;
} HittableList;

typedef struct Sphere {
  const Material *material;
  Point3 center;
  double radius;
} Sphere;

typedef struct MovingSphere {
  const Material *material;
  Point3 center_start, center_end;
  double start, end;
  double radius;
} MovingSphere;

typedef struct BVHNode {
  const Hittable *left;
  const Hittable *right;
  AABB box;
} BVHNode;

struct Hittable {
  HittableType type;
  union {
    HittableList list;
    Sphere sphere;
    MovingSphere moving_sphere;
    BVHNode bvh_node;
  };
};

Hittable *hittable_create_hittable_list(Arena *arena);
Hittable *hittable_create_sphere(Point3 center, double radius,
                                 const Material *material, Arena *arena);
Hittable *hittable_create_moving_sphere(Point3 center_start, Point3 center_end,
                                        double start, double end, double radius,
                                        const Material *material, Arena *arena);
Hittable *hittable_create_bvh_node(const Hittable **objects, size_t start,
                                   size_t end, double time_start,
                                   double time_end, Arena *arena);

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record);
bool hittable_bounding_box(const Hittable *hittable, double time_start,
                           double time_end, AABB *bounding_box);

void hittable_list_add(HittableList *list, const Hittable *hittable,
                       Arena *arena);
bool hittable_list_hit(const HittableList *list, Ray r, double t_min,
                       double t_max, HitRecord *record);
bool hittable_list_bounding_box(const HittableList *list, double time_start,
                                double time_end, AABB *bounding_box);

bool sphere_hit(const Sphere *sphere, Ray r, double t_min, double t_max,
                HitRecord *record);
bool sphere_bounding_box(const Sphere *sphere, double time_start,
                         double time_end, AABB *bounding_box);

Point3 moving_sphere_center(const MovingSphere *sphere, double t);
bool moving_sphere_hit(const MovingSphere *sphere, Ray r, double t_min,
                       double t_max, HitRecord *record);
bool moving_sphere_bounding_box(const MovingSphere *sphere, double time_start,
                                double time_end, AABB *bounding_box);

bool bvh_node_hit(const BVHNode *node, Ray r, double t_min, double t_max,
                  HitRecord *record);
bool bvh_node_bounding_box(const BVHNode *node, double time_start,
                           double time_end, AABB *bounding_box);

#endif /* INCLUDED_HITTABLE_H */

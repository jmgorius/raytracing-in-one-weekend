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
  HITTABLE_XY_RECTANGLE,
  HITTABLE_XZ_RECTANGLE,
  HITTABLE_YZ_RECTANGLE,
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

typedef struct XYRectangle {
  const Material *material;
  double x0, x1, y0, y1, k;
} XYRectangle;

typedef struct XZRectangle {
  const Material *material;
  double x0, x1, z0, z1, k;
} XZRectangle;

typedef struct YZRectangle {
  const Material *material;
  double y0, y1, z0, z1, k;
} YZRectangle;

struct Hittable {
  HittableType type;
  union {
    HittableList list;
    Sphere sphere;
    MovingSphere moving_sphere;
    BVHNode bvh_node;
    XYRectangle xy_rectangle;
    XZRectangle xz_rectangle;
    YZRectangle yz_rectangle;
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
Hittable *hittable_create_xy_rectangle(double x0, double x1, double y0,
                                       double y1, double k,
                                       const Material *material, Arena *arena);
Hittable *hittable_create_xz_rectangle(double x0, double x1, double z0,
                                       double z1, double k,
                                       const Material *material, Arena *arena);
Hittable *hittable_create_yz_rectangle(double y0, double y1, double z0,
                                       double z1, double k,
                                       const Material *material, Arena *arena);

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record);
bool hittable_bounding_box(const Hittable *hittable, double time_start,
                           double time_end, AABB *bounding_box);

void hittable_list_add(HittableList *list, const Hittable *hittable,
                       Arena *arena);

#endif /* INCLUDED_HITTABLE_H */

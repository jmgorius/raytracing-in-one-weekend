#include "hittable.h"
#include "arena.h"
#include "point3.h"
#include "vec3.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void hit_record_set_face_normal(HitRecord *record, Ray r, Vec3 outward_normal) {
  record->front_face = vec3_dot(r.direction, outward_normal) < 0;
  record->normal =
      record->front_face ? outward_normal : vec3_neg(outward_normal);
}

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record) {
  switch (hittable->type) {
  case HITTABLE_LIST:
    return hittable_list_hit((const HittableList *)hittable, r, t_min, t_max,
                             record);
  case HITTABLE_SPHERE:
    return sphere_hit((const Sphere *)hittable, r, t_min, t_max, record);
  case HITTABLE_MOVING_SPHERE:
    return moving_sphere_hit((const MovingSphere *)hittable, r, t_min, t_max,
                             record);
  }
  return false;
}

static void hittable_list_grow(HittableList *list, size_t n, Arena *arena) {
  if (list->objects) {
    const Hittable **new_objects =
        arena_alloc(arena, (list->capacity + n) * sizeof(Hittable *));
    memcpy(new_objects, list->objects, list->size * sizeof(Hittable *));
    list->objects = new_objects;
    list->capacity += n;
  } else {
    list->objects = arena_alloc(arena, n * sizeof(Hittable *));
    list->capacity = n;
    list->size = 0;
  }
}

void hittable_list_add(HittableList *list, const Hittable *hittable,
                       Arena *arena) {
  if (list->capacity == list->size)
    hittable_list_grow(list, list->capacity == 0 ? 16 : list->capacity, arena);
  list->objects[list->size++] = hittable;
}

bool hittable_list_hit(const HittableList *list, Ray r, double t_min,
                       double t_max, HitRecord *record) {
  bool hit_anything = false;
  double closest_so_far = t_max;

  for (size_t i = 0; i < list->size; ++i) {
    if (hittable_hit(list->objects[i], r, t_min, closest_so_far, record)) {
      hit_anything = true;
      closest_so_far = record->t;
    }
  }

  return hit_anything;
}

Sphere *sphere_create(Point3 center, double radius, const Material *material,
                      Arena *arena) {
  Sphere *sphere = arena_alloc(arena, sizeof(Sphere));
  sphere->type = HITTABLE_SPHERE;
  sphere->center = center;
  sphere->radius = radius;
  sphere->material = material;
  return sphere;
}

bool sphere_hit(const Sphere *sphere, Ray r, double t_min, double t_max,
                HitRecord *record) {
  Vec3 oc = point3_sub(r.origin, sphere->center);
  double a = vec3_length2(r.direction);
  double half_b = vec3_dot(oc, r.direction);
  double c = vec3_length2(oc) - sphere->radius * sphere->radius;
  double discriminant = half_b * half_b - a * c;
  if (discriminant < 0)
    return false;

  double square_root = sqrt(discriminant);
  double root = (-half_b - square_root) / a;
  if (root < t_min || t_max < root) {
    root = (-half_b + square_root) / a;
    if (root < t_min || t_max < root)
      return false;
  }

  record->t = root;
  record->p = ray_at(r, root);
  Vec3 outward_normal =
      vec3_div(point3_sub(record->p, sphere->center), sphere->radius);
  hit_record_set_face_normal(record, r, outward_normal);
  record->material = sphere->material;
  return true;
}

MovingSphere *moving_sphere_create(Point3 center_start, Point3 center_end,
                                   double start, double end, double radius,
                                   const Material *material, Arena *arena) {
  assert(start <= end);
  MovingSphere *sphere = arena_alloc(arena, sizeof(MovingSphere));
  sphere->type = HITTABLE_MOVING_SPHERE;
  sphere->center_start = center_start;
  sphere->center_end = center_end;
  sphere->start = start;
  sphere->end = end;
  sphere->radius = radius;
  sphere->material = material;
  return sphere;
}

Point3 moving_sphere_center(const MovingSphere *sphere, double t) {
  Vec3 dir = point3_sub(sphere->center_end, sphere->center_start);
  double c = (t - sphere->start) / (sphere->end - sphere->start);
  return point3_add(sphere->center_start, vec3_mul(c, dir));
}

bool moving_sphere_hit(const MovingSphere *sphere, Ray r, double t_min,
                       double t_max, HitRecord *record) {
  Vec3 oc = point3_sub(r.origin, moving_sphere_center(sphere, r.time));
  double a = vec3_length2(r.direction);
  double half_b = vec3_dot(oc, r.direction);
  double c = vec3_length2(oc) - sphere->radius * sphere->radius;
  double discriminant = half_b * half_b - a * c;
  if (discriminant < 0)
    return false;

  double square_root = sqrt(discriminant);
  double root = (-half_b - square_root) / a;
  if (root < t_min || t_max < root) {
    root = (-half_b + square_root) / a;
    if (root < t_min || t_max < root)
      return false;
  }

  record->t = root;
  record->p = ray_at(r, root);
  Vec3 outward_normal =
      vec3_div(point3_sub(record->p, moving_sphere_center(sphere, r.time)),
               sphere->radius);
  hit_record_set_face_normal(record, r, outward_normal);
  record->material = sphere->material;
  return true;
}

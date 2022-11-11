#include "hittable.h"
#include "point3.h"
#include "vec3.h"

#include <math.h>
#include <stdlib.h>

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
  }
  return false;
}

static void hittable_list_grow(HittableList *list, size_t n) {
  if (list->objects) {
    list->objects = realloc(list->objects, n * sizeof(Hittable *));
    if (!list->objects)
      abort();
    list->capacity = n;
  } else {
    list->objects = malloc(n * sizeof(Hittable *));
    if (!list->objects)
      abort();
    list->capacity = n;
    list->size = 0;
  }
}

void hittable_list_add(HittableList *list, const Hittable *hittable) {
  if (list->capacity == list->size)
    hittable_list_grow(list, list->capacity == 0 ? 16 : list->capacity);
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

void hittable_list_free(HittableList *list) {
  free(list->objects);
  list->objects = 0;
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
  return true;
}

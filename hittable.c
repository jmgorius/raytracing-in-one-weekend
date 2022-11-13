#include "hittable.h"
#include "aabb.h"
#include "arena.h"
#include "point3.h"
#include "utils.h"
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

Hittable *hittable_create_hittable_list(Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_LIST;
  return result;
}

Hittable *hittable_create_sphere(Point3 center, double radius,
                                 const Material *material, Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_SPHERE;
  result->sphere.center = center;
  result->sphere.radius = radius;
  result->sphere.material = material;
  return result;
}

Point3 moving_sphere_center(const MovingSphere *sphere, double t) {
  Vec3 dir = point3_sub(sphere->center_end, sphere->center_start);
  double c = (t - sphere->start) / (sphere->end - sphere->start);
  return point3_add(sphere->center_start, vec3_mul(c, dir));
}

typedef int BoxCompareFunc(const void *lhs, const void *rhs);

#define BOX_COMPARATOR(axis)                                                   \
  static int box_##axis##_compare(const void *lhs, const void *rhs) {          \
    AABB lhs_box, rhs_box;                                                     \
    if (!hittable_bounding_box(*(const Hittable **)lhs, 0, 0, &lhs_box) ||     \
        !hittable_bounding_box(*(const Hittable **)rhs, 0, 0, &rhs_box)) {     \
      fprintf(stderr, "No bounding-box in BVH node");                          \
      exit(1);                                                                 \
    }                                                                          \
                                                                               \
    return lhs_box.min.axis < rhs_box.min.axis;                                \
  }
BOX_COMPARATOR(x)
BOX_COMPARATOR(y)
BOX_COMPARATOR(z)
#undef BOX_COMPARATOR

Hittable *hittable_create_bvh_node(const Hittable **objects, size_t start,
                                   size_t end, double time_start,
                                   double time_end, Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_BVH_NODE;

  int axis = random_int_in_range(0, 2);
  BoxCompareFunc *comparator = (axis == 0)   ? box_x_compare
                               : (axis == 1) ? box_y_compare
                                             : box_z_compare;

  size_t object_span = end - start;
  if (object_span == 1) {
    result->bvh_node.left = result->bvh_node.right = objects[start];
  } else if (object_span == 2) {
    if (comparator(&objects[start], &objects[start + 1])) {
      result->bvh_node.left = objects[start];
      result->bvh_node.right = objects[start + 1];
    } else {
      result->bvh_node.left = objects[start + 1];
      result->bvh_node.right = objects[start];
    }
  } else {
    qsort(objects + start, object_span, sizeof(const Hittable *), comparator);
    size_t mid = start + object_span / 2;
    result->bvh_node.left = hittable_create_bvh_node(
        objects, start, mid, time_start, time_end, arena);
    result->bvh_node.right = hittable_create_bvh_node(
        objects, mid, end, time_start, time_end, arena);
  }

  AABB left_box, right_box;

  if (!hittable_bounding_box(result->bvh_node.left, time_start, time_end,
                             &left_box) ||
      !hittable_bounding_box(result->bvh_node.right, time_start, time_end,
                             &right_box)) {
    fprintf(stderr, "No bounding-box in BVH node");
    exit(1);
  }

  result->bvh_node.box = aabb_surrounding_box(&left_box, &right_box);
  return result;
}

bool hittable_hit(const Hittable *hittable, Ray r, double t_min, double t_max,
                  HitRecord *record) {
  switch (hittable->type) {
  case HITTABLE_LIST:
    return hittable_list_hit(&hittable->list, r, t_min, t_max, record);
  case HITTABLE_SPHERE:
    return sphere_hit(&hittable->sphere, r, t_min, t_max, record);
  case HITTABLE_MOVING_SPHERE:
    return moving_sphere_hit(&hittable->moving_sphere, r, t_min, t_max, record);
  case HITTABLE_BVH_NODE:
    return bvh_node_hit(&hittable->bvh_node, r, t_min, t_max, record);
  }
  return false;
}

bool hittable_bounding_box(const Hittable *hittable, double time_start,
                           double time_end, AABB *bounding_box) {
  switch (hittable->type) {
  case HITTABLE_LIST:
    return hittable_list_bounding_box(&hittable->list, time_start, time_end,
                                      bounding_box);
  case HITTABLE_SPHERE:
    return sphere_bounding_box(&hittable->sphere, time_start, time_end,
                               bounding_box);
  case HITTABLE_MOVING_SPHERE:
    return moving_sphere_bounding_box(&hittable->moving_sphere, time_start,
                                      time_end, bounding_box);
  case HITTABLE_BVH_NODE:
    return bvh_node_bounding_box(&hittable->bvh_node, time_start, time_end,
                                 bounding_box);
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

bool hittable_list_bounding_box(const HittableList *list, double time_start,
                                double time_end, AABB *bounding_box) {
  if (list->size == 0)
    return false;

  AABB temp_box;
  bool first_box = true;

  for (size_t i = 0; i < list->size; ++i) {
    if (!hittable_bounding_box(list->objects[i], time_start, time_end,
                               &temp_box))
      return false;
    *bounding_box =
        first_box ? temp_box : aabb_surrounding_box(bounding_box, &temp_box);
    first_box = false;
  }

  return true;
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

bool sphere_bounding_box(const Sphere *sphere, double time_start,
                         double time_end, AABB *bounding_box) {
  (void)time_start, (void)time_end;
  *bounding_box = (AABB){
      .min = point3_add(sphere->center, (Vec3){-sphere->radius, -sphere->radius,
                                               -sphere->radius}),
      .max = point3_add(sphere->center,
                        (Vec3){sphere->radius, sphere->radius, sphere->radius}),
  };
  return true;
}

Hittable *hittable_create_moving_sphere(Point3 center_start, Point3 center_end,
                                        double start, double end, double radius,
                                        const Material *material,
                                        Arena *arena) {
  assert(start <= end);
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_MOVING_SPHERE;
  result->moving_sphere.center_start = center_start;
  result->moving_sphere.center_end = center_end;
  result->moving_sphere.start = start;
  result->moving_sphere.end = end;
  result->moving_sphere.radius = radius;
  result->moving_sphere.material = material;
  return result;
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

bool moving_sphere_bounding_box(const MovingSphere *sphere, double time_start,
                                double time_end, AABB *bounding_box) {
  AABB box_start = {
      .min =
          point3_add(moving_sphere_center(sphere, time_start),
                     (Vec3){-sphere->radius, -sphere->radius, -sphere->radius}),
      .max = point3_add(moving_sphere_center(sphere, time_start),
                        (Vec3){sphere->radius, sphere->radius, sphere->radius}),
  };
  AABB box_end = {
      .min =
          point3_add(moving_sphere_center(sphere, time_end),
                     (Vec3){-sphere->radius, -sphere->radius, -sphere->radius}),
      .max = point3_add(moving_sphere_center(sphere, time_end),
                        (Vec3){sphere->radius, sphere->radius, sphere->radius}),
  };
  *bounding_box = aabb_surrounding_box(&box_start, &box_end);
  return true;
}

bool bvh_node_hit(const BVHNode *node, Ray r, double t_min, double t_max,
                  HitRecord *record) {
  if (!aabb_hit(&node->box, r, t_min, t_max))
    return false;

  bool hit_left = hittable_hit(node->left, r, t_min, t_max, record);
  bool hit_right =
      hittable_hit(node->right, r, t_min, hit_left ? record->t : t_max, record);

  return hit_left || hit_right;
}

bool bvh_node_bounding_box(const BVHNode *node, double time_start,
                           double time_end, AABB *bounding_box) {
  (void)time_start, (void)time_end;
  *bounding_box = node->box;
  return true;
}

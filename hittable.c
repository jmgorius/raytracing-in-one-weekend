#include "hittable.h"
#include "aabb.h"
#include "arena.h"
#include "material.h"
#include "point3.h"
#include "ray.h"
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

static Point3 moving_sphere_center(const MovingSphere *sphere, double t) {
  Vec3 dir = point3_sub(sphere->center_end, sphere->center_start);
  double c = (t - sphere->start) / (sphere->end - sphere->start);
  return point3_add(sphere->center_start, vec3_mul(c, dir));
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

Hittable *hittable_create_xy_rectangle(double x0, double x1, double y0,
                                       double y1, double k,
                                       const Material *material, Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_XY_RECTANGLE;
  result->xy_rectangle.x0 = x0;
  result->xy_rectangle.x1 = x1;
  result->xy_rectangle.y0 = y0;
  result->xy_rectangle.y1 = y1;
  result->xy_rectangle.k = k;
  result->xy_rectangle.material = material;
  return result;
}

Hittable *hittable_create_xz_rectangle(double x0, double x1, double z0,
                                       double z1, double k,
                                       const Material *material, Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_XZ_RECTANGLE;
  result->xz_rectangle.x0 = x0;
  result->xz_rectangle.x1 = x1;
  result->xz_rectangle.z0 = z0;
  result->xz_rectangle.z1 = z1;
  result->xz_rectangle.k = k;
  result->xz_rectangle.material = material;
  return result;
}

Hittable *hittable_create_yz_rectangle(double y0, double y1, double z0,
                                       double z1, double k,
                                       const Material *material, Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_YZ_RECTANGLE;
  result->yz_rectangle.y0 = y0;
  result->yz_rectangle.y1 = y1;
  result->yz_rectangle.z0 = z0;
  result->yz_rectangle.z1 = z1;
  result->yz_rectangle.k = k;
  result->yz_rectangle.material = material;
  return result;
}

Hittable *hittable_create_box(Point3 p0, Point3 p1, const Material *material,
                              Arena *arena) {
  Hittable *result = arena_alloc(arena, sizeof(Hittable));
  result->type = HITTABLE_BOX;
  result->box.min = p0;
  result->box.max = p1;
  hittable_list_add(&result->box.sides,
                    hittable_create_xy_rectangle(p0.x, p1.x, p0.y, p1.y, p0.z,
                                                 material, arena),
                    arena);
  hittable_list_add(&result->box.sides,
                    hittable_create_xy_rectangle(p0.x, p1.x, p0.y, p1.y, p1.z,
                                                 material, arena),
                    arena);
  hittable_list_add(&result->box.sides,
                    hittable_create_xz_rectangle(p0.x, p1.x, p0.z, p1.z, p0.y,
                                                 material, arena),
                    arena);
  hittable_list_add(&result->box.sides,
                    hittable_create_xz_rectangle(p0.x, p1.x, p0.z, p1.z, p1.y,
                                                 material, arena),
                    arena);
  hittable_list_add(&result->box.sides,
                    hittable_create_yz_rectangle(p0.y, p1.y, p0.z, p1.z, p0.x,
                                                 material, arena),
                    arena);
  hittable_list_add(&result->box.sides,
                    hittable_create_yz_rectangle(p0.y, p1.y, p0.z, p1.z, p1.x,
                                                 material, arena),
                    arena);
  return result;
}

static bool hittable_list_hit(const HittableList *list, Ray r, double t_min,
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

static void get_sphere_uv(Point3 p, double *u, double *v) {
  double theta = acos(-p.y);
  double phi = atan2(-p.z, p.x) + M_PI;
  *u = phi / (2 * M_PI);
  *v = theta / M_PI;
}

static bool sphere_hit(const Sphere *sphere, Ray r, double t_min, double t_max,
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
  get_sphere_uv((Point3){outward_normal.x, outward_normal.y, outward_normal.z},
                &record->u, &record->v);
  record->material = sphere->material;
  return true;
}

static bool moving_sphere_hit(const MovingSphere *sphere, Ray r, double t_min,
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
  get_sphere_uv((Point3){outward_normal.x, outward_normal.y, outward_normal.z},
                &record->u, &record->v);
  record->material = sphere->material;
  return true;
}

static bool bvh_node_hit(const BVHNode *node, Ray r, double t_min, double t_max,
                         HitRecord *record) {
  if (!aabb_hit(&node->box, r, t_min, t_max))
    return false;

  bool hit_left = hittable_hit(node->left, r, t_min, t_max, record);
  bool hit_right =
      hittable_hit(node->right, r, t_min, hit_left ? record->t : t_max, record);

  return hit_left || hit_right;
}

static bool xy_rectangle_hit(const XYRectangle *rectangle, Ray r, double t_min,
                             double t_max, HitRecord *record) {
  double t = (rectangle->k - r.origin.z) / r.direction.z;
  if (t < t_min || t > t_max)
    return false;

  double x = r.origin.x + t * r.direction.x;
  double y = r.origin.y + t * r.direction.y;

  if (x < rectangle->x0 || x > rectangle->x1 || y < rectangle->y0 ||
      y > rectangle->y1)
    return false;

  record->u = (x - rectangle->x0) / (rectangle->x1 - rectangle->x0);
  record->v = (y - rectangle->y0) / (rectangle->y1 - rectangle->y0);
  record->t = t;

  Vec3 outward_normal = {0.0, 0.0, 1.0};
  hit_record_set_face_normal(record, r, outward_normal);
  record->material = rectangle->material;
  record->p = ray_at(r, t);

  return true;
}

static bool xz_rectangle_hit(const XZRectangle *rectangle, Ray r, double t_min,
                             double t_max, HitRecord *record) {
  double t = (rectangle->k - r.origin.y) / r.direction.y;
  if (t < t_min || t > t_max)
    return false;

  double x = r.origin.x + t * r.direction.x;
  double z = r.origin.z + t * r.direction.z;

  if (x < rectangle->x0 || x > rectangle->x1 || z < rectangle->z0 ||
      z > rectangle->z1)
    return false;

  record->u = (x - rectangle->x0) / (rectangle->x1 - rectangle->x0);
  record->v = (z - rectangle->z0) / (rectangle->z1 - rectangle->z0);
  record->t = t;

  Vec3 outward_normal = {0.0, 1.0, 0.0};
  hit_record_set_face_normal(record, r, outward_normal);
  record->material = rectangle->material;
  record->p = ray_at(r, t);

  return true;
}

static bool yz_rectangle_hit(const YZRectangle *rectangle, Ray r, double t_min,
                             double t_max, HitRecord *record) {
  double t = (rectangle->k - r.origin.x) / r.direction.x;
  if (t < t_min || t > t_max)
    return false;

  double y = r.origin.y + t * r.direction.y;
  double z = r.origin.z + t * r.direction.z;

  if (y < rectangle->y0 || y > rectangle->y1 || z < rectangle->z0 ||
      z > rectangle->z1)
    return false;

  record->u = (y - rectangle->y0) / (rectangle->y1 - rectangle->y0);
  record->v = (z - rectangle->z0) / (rectangle->z1 - rectangle->z0);
  record->t = t;

  Vec3 outward_normal = {1.0, 0.0, 0.0};
  hit_record_set_face_normal(record, r, outward_normal);
  record->material = rectangle->material;
  record->p = ray_at(r, t);

  return true;
}

static bool box_hit(const Box *box, Ray r, double t_min, double t_max,
                    HitRecord *record) {
  return hittable_list_hit(&box->sides, r, t_min, t_max, record);
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
  case HITTABLE_XY_RECTANGLE:
    return xy_rectangle_hit(&hittable->xy_rectangle, r, t_min, t_max, record);
  case HITTABLE_XZ_RECTANGLE:
    return xz_rectangle_hit(&hittable->xz_rectangle, r, t_min, t_max, record);
  case HITTABLE_YZ_RECTANGLE:
    return yz_rectangle_hit(&hittable->yz_rectangle, r, t_min, t_max, record);
  case HITTABLE_BOX:
    return box_hit(&hittable->box, r, t_min, t_max, record);
  }
  return false;
}

static bool hittable_list_bounding_box(const HittableList *list,
                                       double time_start, double time_end,
                                       AABB *bounding_box) {
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

static bool sphere_bounding_box(const Sphere *sphere, AABB *bounding_box) {
  *bounding_box = (AABB){
      .min = point3_add(sphere->center, (Vec3){-sphere->radius, -sphere->radius,
                                               -sphere->radius}),
      .max = point3_add(sphere->center,
                        (Vec3){sphere->radius, sphere->radius, sphere->radius}),
  };
  return true;
}

static bool moving_sphere_bounding_box(const MovingSphere *sphere,
                                       double time_start, double time_end,
                                       AABB *bounding_box) {
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

static bool bvh_node_bounding_box(const BVHNode *node, AABB *bounding_box) {
  *bounding_box = node->box;
  return true;
}

static bool xy_rectangle_bounding_box(const XYRectangle *rectangle,
                                      AABB *bounding_box) {
  /* Pad the bounding box to make sure it is not zero-width */
  *bounding_box = (AABB){
      .min = {rectangle->x0, rectangle->y0, rectangle->k - 0.0001},
      .max = {rectangle->x1, rectangle->y1, rectangle->k + 0.0001},
  };
  return true;
}

static bool xz_rectangle_bounding_box(const XZRectangle *rectangle,
                                      AABB *bounding_box) {
  /* Pad the bounding box to make sure it is not zero-width */
  *bounding_box = (AABB){
      .min = {rectangle->x0, rectangle->k - 0.0001, rectangle->z0},
      .max = {rectangle->x1, rectangle->k + 0.0001, rectangle->z1},
  };
  return true;
}

static bool yz_rectangle_bounding_box(const YZRectangle *rectangle,
                                      AABB *bounding_box) {
  /* Pad the bounding box to make sure it is not zero-width */
  *bounding_box = (AABB){
      .min = {rectangle->k - 0.0001, rectangle->y0, rectangle->z0},
      .max = {rectangle->k + 0.0001, rectangle->y1, rectangle->z1},
  };
  return true;
}

static bool box_bounding_box(const Box *box, AABB *bounding_box) {
  *bounding_box = (AABB){box->min, box->max};
  return true;
}

bool hittable_bounding_box(const Hittable *hittable, double time_start,
                           double time_end, AABB *bounding_box) {
  switch (hittable->type) {
  case HITTABLE_LIST:
    return hittable_list_bounding_box(&hittable->list, time_start, time_end,
                                      bounding_box);
  case HITTABLE_SPHERE:
    return sphere_bounding_box(&hittable->sphere, bounding_box);
  case HITTABLE_MOVING_SPHERE:
    return moving_sphere_bounding_box(&hittable->moving_sphere, time_start,
                                      time_end, bounding_box);
  case HITTABLE_BVH_NODE:
    return bvh_node_bounding_box(&hittable->bvh_node, bounding_box);
  case HITTABLE_XY_RECTANGLE:
    return xy_rectangle_bounding_box(&hittable->xy_rectangle, bounding_box);
  case HITTABLE_XZ_RECTANGLE:
    return xz_rectangle_bounding_box(&hittable->xz_rectangle, bounding_box);
  case HITTABLE_YZ_RECTANGLE:
    return yz_rectangle_bounding_box(&hittable->yz_rectangle, bounding_box);
  case HITTABLE_BOX:
    return box_bounding_box(&hittable->box, bounding_box);
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

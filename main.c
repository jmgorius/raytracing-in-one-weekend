#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "arena.h"
#include "camera.h"
#include "color.h"
#include "hittable.h"
#include "material.h"
#include "point3.h"
#include "ray.h"
#include "utils.h"
#include "vec3.h"

static Color ray_color(Ray r, const Hittable *world, int depth) {
  if (depth <= 0)
    return (Color){0, 0, 0};

  HitRecord record;
  if (hittable_hit(world, r, 0.001, DBL_MAX, &record)) {
    Ray scattered;
    Color attenuation;
    if (material_scatter(record.material, r, &record, &attenuation, &scattered))
      return color_mul(attenuation, ray_color(scattered, world, depth - 1));
    return (Color){0, 0, 0};
  }

  Vec3 unit_direction = vec3_normalize(r.direction);
  double t = 0.5 * (unit_direction.y + 1.0);
  Color gradient1 = {1.0, 1.0, 1.0};
  Color gradient2 = {0.5, 0.7, 1.0};
  return color_lerp(gradient1, gradient2, t);
}

static const Hittable *generate_random_scene(Arena *arena) {
  static HittableList world = {.type = HITTABLE_LIST};

  static Lambertian ground_material = {.type = MATERIAL_LAMBERTIAN,
                                       .albedo = (Color){0.5, 0.5, 0.5}};
  Sphere *ground_sphere = arena_alloc(arena, sizeof(Sphere));
  ground_sphere->type = HITTABLE_SPHERE;
  ground_sphere->center = (Point3){0.0, -1000.0, 0.0};
  ground_sphere->radius = 1000.0;
  ground_sphere->material = (const Material *)&ground_material;
  hittable_list_add(&world, (const Hittable *)ground_sphere, arena);

  static Lambertian lambertian = {.type = MATERIAL_LAMBERTIAN,
                                  .albedo = {0.4, 0.2, 0.1}};
  static Metal metal = {
      .type = MATERIAL_METAL, .albedo = {0.7, 0.6, 0.5}, .fuzziness = 0.0};
  static Dielectric glass = {.type = MATERIAL_DIELECTRIC, .eta = 1.5};

  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      double choose_material = random_double();
      Point3 center = {a + 0.9 * random_double(), 0.2,
                       b + 0.9 * random_double()};

      if (vec3_length(point3_sub(center, (Point3){4.0, 0.2, 0.0})) > 0.9) {
        if (choose_material < 0.8) {
          Color albedo = color_mul(color_random(), color_random());
          Lambertian *material = arena_alloc(arena, sizeof(Lambertian));
          material->type = MATERIAL_LAMBERTIAN;
          material->albedo = albedo;

          Sphere *sphere = arena_alloc(arena, sizeof(Sphere));
          sphere->type = HITTABLE_SPHERE;
          sphere->center = center;
          sphere->radius = 0.2;
          sphere->material = (const Material *)material;
          hittable_list_add(&world, (const Hittable *)sphere, arena);
        } else if (choose_material < 0.95) {
          Color albedo = color_random_in_range(0.5, 1);
          double fuzziness = random_double_in_range(0.5, 1.0);
          Metal *material = arena_alloc(arena, sizeof(Metal));
          material->type = MATERIAL_METAL;
          material->albedo = albedo;
          material->fuzziness = fuzziness;

          Sphere *sphere = arena_alloc(arena, sizeof(Sphere));
          sphere->type = HITTABLE_SPHERE;
          sphere->center = center;
          sphere->radius = 0.2;
          sphere->material = (const Material *)material;
          hittable_list_add(&world, (const Hittable *)sphere, arena);
        } else {
          Sphere *sphere = arena_alloc(arena, sizeof(Sphere));
          sphere->type = HITTABLE_SPHERE;
          sphere->center = center;
          sphere->radius = 0.2;
          sphere->material = (const Material *)&glass;
          hittable_list_add(&world, (const Hittable *)sphere, arena);
        }
      }
    }
  }

  Sphere *sphere1 = arena_alloc(arena, sizeof(Sphere));
  sphere1->type = HITTABLE_SPHERE;
  sphere1->center = (Point3){0.0, 1.0, 0.0};
  sphere1->radius = 1.0;
  sphere1->material = (const Material *)&glass;
  hittable_list_add(&world, (const Hittable *)sphere1, arena);
  Sphere *sphere2 = arena_alloc(arena, sizeof(Sphere));
  sphere2->type = HITTABLE_SPHERE;
  sphere2->center = (Point3){-4.0, 1.0, 0.0};
  sphere2->radius = 1.0;
  sphere2->material = (const Material *)&lambertian;
  hittable_list_add(&world, (const Hittable *)sphere2, arena);
  Sphere *sphere3 = arena_alloc(arena, sizeof(Sphere));
  sphere3->type = HITTABLE_SPHERE;
  sphere3->center = (Point3){4.0, 1.0, 0.0};
  sphere3->radius = 1.0;
  sphere3->material = (const Material *)&metal;
  hittable_list_add(&world, (const Hittable *)sphere3, arena);

  return (const Hittable *)&world;
}

int main(void) {
  /* Memory management */

  const unsigned int buffer_size = 1 * 1024 * 1024;
  void *buffer = malloc(buffer_size);
  if (!buffer)
    abort();
  Arena arena = {0};
  arena_init(&arena, buffer, buffer_size);

  /* Image parameters */

  const double aspect_ratio = 3.0 / 2.0;
  const int image_width = 1200;
  const int image_height = (int)(image_width / aspect_ratio);
  const int samples_per_pixel = 500;
  const int max_depth = 50;

  /* World */

  const Hittable *world = generate_random_scene(&arena);

  /* Camera */

  Point3 look_from = {13.0, 2.0, 3.0};
  Point3 look_at = {0.0, 0.0, 0.0};
  Vec3 up = {0.0, 1.0, 0.0};
  double dist_to_focus = 10.0;
  double aperture = 0.1;

  Camera camera;
  camera_init(&camera, look_from, look_at, up, 20, aspect_ratio, aperture,
              dist_to_focus);

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      Color pixel_color = {0, 0, 0};
      for (int s = 0; s < samples_per_pixel; ++s) {
        double u = (i + random_double()) / (image_width - 1);
        double v = (j + random_double()) / (image_height - 1);
        Ray r = camera_get_ray(&camera, u, v);
        pixel_color = color_add(pixel_color, ray_color(r, world, max_depth));
      }
      color_write(stdout, pixel_color, samples_per_pixel);
    }
  }

  fprintf(stderr, "\nDone.\n");

  free(buffer);

  return 0;
}

#include "arena.c"
#include "camera.c"
#include "color.c"
#include "hittable.c"
#include "material.c"
#include "point3.c"
#include "ray.c"
#include "utils.c"
#include "vec3.c"

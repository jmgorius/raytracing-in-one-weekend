#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "camera.h"
#include "color.h"
#include "hittable.h"
#include "material.h"
#include "point3.h"
#include "ray.h"
#include "utils.h"
#include "vec3.h"

Color ray_color(Ray r, const Hittable *world, int depth) {
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

int main(void) {
  /* Image parameters */
  const double aspect_ratio = 16.0 / 9.0;
  const int image_width = 256;
  const int image_height = (int)(image_width / aspect_ratio);
  const int samples_per_pixel = 100;
  const int max_depth = 50;

  /* World */

  HittableList world = {.type = HITTABLE_LIST};

  Lambertian material_ground = {.type = MATERIAL_LAMBERTIAN,
                                .albedo = (Color){0.8, 0.8, 0.0}};
  Lambertian material_center = {.type = MATERIAL_LAMBERTIAN,
                                .albedo = (Color){0.7, 0.3, 0.3}};
  Metal material_left = {.type = MATERIAL_METAL,
                         .albedo = (Color){0.8, 0.8, 0.8},
                         .fuzziness = 0.3};
  Metal material_right = {.type = MATERIAL_METAL,
                          .albedo = (Color){0.8, 0.6, 0.2},
                          .fuzziness = 1.0};

  Sphere sphere_ground = {
      .type = HITTABLE_SPHERE,
      .center = (Point3){0.0, -100.5, -1},
      .radius = 100.0,
      .material = (const Material *)&material_ground,
  };
  Sphere sphere_center = {
      .type = HITTABLE_SPHERE,
      .center = (Point3){0.0, 0.0, -1.0},
      .radius = 0.5,
      .material = (const Material *)&material_center,
  };
  Sphere sphere_left = {
      .type = HITTABLE_SPHERE,
      .center = (Point3){-1.0, 0.0, -1.0},
      .radius = 0.5,
      .material = (const Material *)&material_left,
  };
  Sphere sphere_right = {
      .type = HITTABLE_SPHERE,
      .center = (Point3){1.0, 0.0, -1.0},
      .radius = 0.5,
      .material = (const Material *)&material_right,
  };
  hittable_list_add(&world, (const Hittable *)&sphere_ground);
  hittable_list_add(&world, (const Hittable *)&sphere_center);
  hittable_list_add(&world, (const Hittable *)&sphere_left);
  hittable_list_add(&world, (const Hittable *)&sphere_right);

  /* Camera */
  Camera camera;
  camera_init(&camera, aspect_ratio);

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      Color pixel_color = {0, 0, 0};
      for (int s = 0; s < samples_per_pixel; ++s) {
        double u = (i + random_double()) / (image_width - 1);
        double v = (j + random_double()) / (image_height - 1);
        Ray r = camera_get_ray(&camera, u, v);
        pixel_color = color_add(
            pixel_color, ray_color(r, (const Hittable *)&world, max_depth));
      }
      color_write(stdout, pixel_color, samples_per_pixel);
    }
  }

  fprintf(stderr, "\nDone.\n");

  hittable_list_free(&world);

  return 0;
}

#include "camera.c"
#include "color.c"
#include "hittable.c"
#include "material.c"
#include "point3.c"
#include "ray.c"
#include "utils.c"
#include "vec3.c"

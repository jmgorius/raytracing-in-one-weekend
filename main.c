#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "arena.h"
#include "camera.h"
#include "color.h"
#include "hittable.h"
#include "material.h"
#include "point3.h"
#include "ray.h"
#include "scenes.h"
#include "utils.h"
#include "vec3.h"

#ifndef SCENE_SELECT
#define SCENE_SELECT SCENE_COMPLEX
#endif

static Color ray_color(Ray r, Color background_color, const Hittable *world,
                       int depth) {
  if (depth <= 0)
    return (Color){0, 0, 0};

  HitRecord record;
  if (!hittable_hit(world, r, 0.001, DBL_MAX, &record))
    return background_color;

  Ray scattered;
  Color attenuation;
  Color emitted =
      material_emitted(record.material, record.u, record.v, record.p);

  if (!material_scatter(record.material, r, &record, &attenuation, &scattered))
    return emitted;

  return color_add(emitted,
                   color_mul(attenuation, ray_color(scattered, background_color,
                                                    world, depth - 1)));
}

int main(int argc, char *argv[]) {
  srand(time(0));

  if (argc < 2) {
    fprintf(stderr, "Expected output file name\n");
    return 1;
  }

  FILE *output_file = fopen(argv[1], "wb");
  if (!output_file) {
    fprintf(stderr, "Failed to open output file %s: %s", argv[1],
            strerror(errno));
    return 1;
  }

  /* Memory management */

  const unsigned int buffer_size = 1 * 1024 * 1024;
  void *buffer = malloc(buffer_size);
  if (!buffer)
    abort();
  Arena arena = {0};
  arena_init(&arena, buffer, buffer_size);

  /* Image parameters */

  double aspect_ratio = 16.0 / 9.0;
  int image_width = 400;
  int samples_per_pixel = 400;
  const int max_depth = 50;

  Point3 look_from = {0.0, 0.0, 1.0};
  Point3 look_at = {0.0, 0.0, 0.0};
  double vfov = 40.0;
  double aperture = 0.0;
  double dist_to_focus = 10.0;
  Color background_color = {0.0, 0.0, 0.0};

  /* World */

  const Hittable *world = 0;

#if SCENE_SELECT == SCENE_RANDOM
  world = random_scene(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13, 2, 3};
  look_at = (Point3){0, 0, 0};
  vfov = 20.0;
  aperture = 0.1;
#elif SCENE_SELECT == SCENE_TWO_SPHERES
  world = two_spheres(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13, 2, 3};
  look_at = (Point3){0, 0, 0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_TWO_PERLIN_SPHERES
  world = two_perlin_spheres(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13, 2, 3};
  look_at = (Point3){0, 0, 0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_EARTH
  world = earth(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13, 2, 3};
  look_at = (Point3){0, 0, 0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_SIMPLE_LIGHT
  world = simple_light(&arena);
  background_color = (Color){0.0, 0.0, 0.0};
  look_from = (Point3){26, 3, 6};
  look_at = (Point3){0, 2, 0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_CORNELL_BOX
  world = cornell_box(&arena);
  aspect_ratio = 1.0;
  image_width = 600;
  samples_per_pixel = 1000;
  background_color = (Color){0.0, 0.0, 0.0};
  look_from = (Point3){278, 278, -800};
  look_at = (Point3){278, 278, 0};
  vfov = 40.0;
#elif SCENE_SELECT == SCENE_CORNELL_BOX_SMOKE
  world = cornell_box_smoke(&arena);
  aspect_ratio = 1.0;
  image_width = 600;
  samples_per_pixel = 200;
  look_from = (Point3){278, 278, -800};
  look_at = (Point3){278, 278, 0};
  vfov = 40.0;
#elif SCENE_SELECT == SCENE_COMPLEX
  world = complex_scene(&arena);
  aspect_ratio = 1.0;
  image_width = 800;
  samples_per_pixel = 10000;
  background_color = (Color){0.0, 0.0, 0.0};
  look_from = (Point3){478, 278, -600};
  look_at = (Point3){278, 278, 0};
  vfov = 40.0;
#else
#error Unknown scene selected
#endif

  int image_height = (int)(image_width / aspect_ratio);
  Vec3 up = {0.0, 1.0, 0.0};

  Camera camera;
  camera_init(&camera, look_from, look_at, up, vfov, aspect_ratio, aperture,
              dist_to_focus, 0.0, 1.0);

  fprintf(output_file, "P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      Color pixel_color = {0, 0, 0};
      for (int s = 0; s < samples_per_pixel; ++s) {
        double u = (i + random_double()) / (image_width - 1);
        double v = (j + random_double()) / (image_height - 1);
        Ray r = camera_get_ray(&camera, u, v);
        pixel_color = color_add(
            pixel_color, ray_color(r, background_color, world, max_depth));
      }
      color_write(output_file, pixel_color, samples_per_pixel);
    }
  }

  fprintf(stderr, "\nDone.\n");

  free(buffer);
  fclose(output_file);

  return 0;
}

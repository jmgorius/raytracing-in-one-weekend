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
#include "perlin.h"
#include "point3.h"
#include "ray.h"
#include "texture.h"
#include "utils.h"
#include "vec3.h"

#define SCENE_RANDOM 0
#define SCENE_TWO_SPHERES 1
#define SCENE_TWO_PERLIN_SPHERES 2
#define SCENE_EARTH 3
#define SCENE_SIMPLE_LIGHT 4

#ifndef SCENE_SELECT
#define SCENE_SELECT SCENE_SIMPLE_LIGHT
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

static const Hittable *generate_random_scene(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  const Texture *checker = texture_create_checker_solid_color(
      (Color){0.2, 0.3, 0.1}, (Color){0.9, 0.9, 0.9}, 10.0, arena);
  const Material *ground_material = material_create_lambertian(checker, arena);
  const Hittable *ground_sphere = hittable_create_sphere(
      (Point3){0.0, -1000.0, 0.0}, 1000.0, ground_material, arena);
  hittable_list_add(&world->list, ground_sphere, arena);

  const Material *glass = material_create_dielectric(1.5, arena);

  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      double choose_material = random_double();
      Point3 center = {a + 0.9 * random_double(), 0.2,
                       b + 0.9 * random_double()};

      if (vec3_length(point3_sub(center, (Point3){4.0, 0.2, 0.0})) > 0.9) {
        if (choose_material < 0.8) {
          const Material *material = material_create_lambertian_color(
              color_mul(color_random(), color_random()), arena);
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, material, arena);
          hittable_list_add(&world->list, sphere, arena);
        } else if (choose_material < 0.95) {
          const Material *material = material_create_metal_color(
              color_random_in_range(0.5, 1), random_double_in_range(0.5, 1.0),
              arena);
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, material, arena);
          hittable_list_add(&world->list, sphere, arena);
        } else {
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, glass, arena);
          hittable_list_add(&world->list, sphere, arena);
        }
      }
    }
  }

  const Material *lambertian =
      material_create_lambertian_color((Color){0.4, 0.2, 0.1}, arena);
  const Material *metal =
      material_create_metal_color((Color){0.7, 0.6, 0.5}, 0.0, arena);

  const Hittable *sphere1 =
      hittable_create_sphere((Point3){0.0, 1.0, 0.0}, 1.0, glass, arena);
  hittable_list_add(&world->list, sphere1, arena);
  const Hittable *sphere2 =
      hittable_create_sphere((Point3){-4.0, 1.0, 0.0}, 1.0, lambertian, arena);
  hittable_list_add(&world->list, sphere2, arena);
  const Hittable *sphere3 =
      hittable_create_sphere((Point3){4.0, 1.0, 0.0}, 1.0, metal, arena);
  hittable_list_add(&world->list, sphere3, arena);

  Hittable *bvh_root = hittable_create_bvh_node(
      world->list.objects, 0, world->list.size, 0.0, 1.0, arena);
  return bvh_root;
}

static Hittable *two_spheres(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *checker = texture_create_checker_solid_color(
      (Color){0.2, 0.3, 0.1}, (Color){0.9, 0.9, 0.9}, 10.0, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -10.0, 0.0}, 10.0,
                             material_create_lambertian(checker, arena), arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 10.0, 0.0}, 10.0,
                             material_create_lambertian(checker, arena), arena),
      arena);

  return world;
}

static Hittable *two_perlin_spheres(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *perlin_texture =
      texture_create_perlin_noise(4.0, PERLIN_DEFAULT_POINT_COUNT, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -1000.0, 0.0}, 1000.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 2.0, 0.0}, 2.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);

  return world;
}

static Hittable *earth(Arena *arena) {
  Texture *earth_texture = texture_create_image("assets/earthmap.jpg", arena);
  Material *earth_surface = material_create_lambertian(earth_texture, arena);
  Hittable *globe = hittable_create_sphere((Point3){0.0, 0.0, 0.0}, 2.0,
                                           earth_surface, arena);
  return globe;
}

static Hittable *simple_light(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *perlin_texture =
      texture_create_perlin_noise(4.0, PERLIN_DEFAULT_POINT_COUNT, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -1000.0, 0.0}, 1000.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 2.0, 0.0}, 2.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);

  Material *diffuse_light =
      material_create_diffuse_light_color((Color){4.0, 4.0, 4.0}, arena);
  hittable_list_add(&world->list,
                    hittable_create_xy_rectangle(3.0, 5.0, 1.0, 3.0, -2.0,
                                                 diffuse_light, arena),
                    arena);

  return world;
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

  const double aspect_ratio = 16.0 / 9.0;
  const int image_width = 400;
  const int image_height = (int)(image_width / aspect_ratio);
  const int samples_per_pixel = 400;
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
  world = generate_random_scene(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13.0, 2.0, 3.0};
  look_at = (Point3){0.0, 0.0, 0.0};
  vfov = 20.0;
  aperture = 0.1;
#elif SCENE_SELECT == SCENE_TWO_SPHERES
  world = two_spheres(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13.0, 2.0, 3.0};
  look_at = (Point3){0.0, 0.0, 0.0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_TWO_PERLIN_SPHERES
  world = two_perlin_spheres(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13.0, 2.0, 3.0};
  look_at = (Point3){0.0, 0.0, 0.0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_EARTH
  world = earth(&arena);
  background_color = (Color){0.7, 0.8, 1.0};
  look_from = (Point3){13.0, 2.0, 3.0};
  look_at = (Point3){0.0, 0.0, 0.0};
  vfov = 20.0;
#elif SCENE_SELECT == SCENE_SIMPLE_LIGHT
  world = simple_light(&arena);
  background_color = (Color){0.0, 0.0, 0.0};
  look_from = (Point3){26.0, 3.0, 6.0};
  look_at = (Point3){0.0, 2.0, 0.0};
  vfov = 20.0;
#endif

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

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "camera.h"
#include "color.h"
#include "hittable.h"
#include "point3.h"
#include "ray.h"
#include "utils.h"
#include "vec3.h"

Color ray_color(Ray r, Hittable world, int depth) {
  if (depth <= 0)
    return (Color){0, 0, 0};

  HitRecord record;
  if (hittable_hit(&world, r, 0.001, DBL_MAX, &record)) {
    Point3 target = point3_add(
        record.p, vec3_add(record.normal, vec3_random_in_unit_sphere()));
    return color_mul(0.5,
                     ray_color((Ray){record.p, point3_sub(target, record.p)},
                               world, depth - 1));
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
  HittableList object_list = {0};
  hittable_list_add(&object_list, make_hittable_sphere(&(Sphere){
                                      .center = (Point3){0, 0, -1},
                                      .radius = 0.5,
                                  }));
  hittable_list_add(&object_list, make_hittable_sphere(&(Sphere){
                                      .center = (Point3){0, -100.5, -1},
                                      .radius = 100,
                                  }));
  Hittable world = make_hittable_list(&object_list);

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
        pixel_color = color_add(pixel_color, ray_color(r, world, max_depth));
      }
      color_write(stdout, pixel_color, samples_per_pixel);
    }
  }

  fprintf(stderr, "\nDone.\n");

  hittable_list_free(&object_list);

  return 0;
}

#include "camera.c"
#include "color.c"
#include "hittable.c"
#include "point3.c"
#include "ray.c"
#include "utils.c"
#include "vec3.c"

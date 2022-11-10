#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "color.h"
#include "hittable.h"
#include "point3.h"
#include "ray.h"
#include "vec3.h"

Color ray_color(Ray r, Hittable world) {
  HitRecord record;
  if (hittable_hit(&world, r, 0, DBL_MAX, &record)) {
    return (Color){0.5 * (record.normal.x + 1), 0.5 * (record.normal.y + 1),
                   0.5 * (record.normal.z + 1)};
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

  /* World */
  HittableList world = {0};
  hittable_list_add(&world, make_hittable_sphere(&(Sphere){
                                .center = (Point3){0, 0, -1},
                                .radius = 0.5,
                            }));
  hittable_list_add(&world, make_hittable_sphere(&(Sphere){
                                .center = (Point3){0, -100.5, -1},
                                .radius = 100,
                            }));

  /* Camera parameters */
  double viewport_height = 2;
  double viewport_width = aspect_ratio * viewport_height;
  double focal_length = 1.0;

  Point3 origin = {0};
  Vec3 horizontal = {viewport_width, 0, 0};
  Vec3 vertical = {0, viewport_height, 0};
  Vec3 offset = vec3_add(vec3_div(horizontal, 2), vec3_div(vertical, 2));
  offset = vec3_add(offset, (Vec3){0, 0, focal_length});
  Point3 lower_left_corner = point3_add(origin, vec3_neg(offset));

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      double u = (double)i / (image_width - 1);
      double v = (double)j / (image_height - 1);
      Point3 screen_point =
          point3_add(lower_left_corner,
                     vec3_add(vec3_mul(u, horizontal), vec3_mul(v, vertical)));
      Vec3 direction = point3_sub(screen_point, origin);
      Ray r = {origin, direction};
      Color pixel_color = ray_color(r, make_hittable_list(&world));
      color_write(stdout, pixel_color);
    }
  }

  fprintf(stderr, "\nDone.\n");

  hittable_list_free(&world);

  return 0;
}

#include "color.c"
#include "hittable.c"
#include "point3.c"
#include "ray.c"
#include "utils.c"
#include "vec3.c"

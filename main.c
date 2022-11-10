#include <stdio.h>

#include "color.h"
#include "point3.h"
#include "ray.h"
#include "vec3.h"

color ray_color(ray r) {
  vec3 unit_direction = vec3_normalize(r.direction);
  double t = 0.5 * (unit_direction.y + 1.0);
  color gradient1 = {1.0, 1.0, 1.0};
  color gradient2 = {0.5, 0.7, 1.0};
  return color_lerp(gradient1, gradient2, t);
}

int main(void) {
  /* Image parameters */
  const double aspect_ratio = 16.0 / 9.0;
  const int image_width = 256;
  const int image_height = (int)(image_width / aspect_ratio);

  /* Camera parameters */
  double viewport_height = 2;
  double viewport_width = aspect_ratio * viewport_height;
  double focal_length = 1.0;

  point3 origin = {0};
  vec3 horizontal = {viewport_width, 0, 0};
  vec3 vertical = {0, viewport_height, 0};
  vec3 offset = vec3_add(vec3_div(horizontal, 2), vec3_div(vertical, 2));
  offset = vec3_add(offset, (vec3){0, 0, focal_length});
  point3 lower_left_corner = point3_add(origin, vec3_neg(offset));

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      double u = (double)i / (image_width - 1);
      double v = (double)j / (image_height - 1);
      point3 screen_point =
          point3_add(lower_left_corner,
                     vec3_add(vec3_mul(u, horizontal), vec3_mul(v, vertical)));
      vec3 direction = point3_sub(screen_point, origin);
      ray r = {origin, direction};
      color pixel_color = ray_color(r);
      color_write(stdout, pixel_color);
    }
  }

  fprintf(stderr, "\nDone.\n");

  return 0;
}

#include "color.c"
#include "point3.c"
#include "ray.c"
#include "vec3.c"

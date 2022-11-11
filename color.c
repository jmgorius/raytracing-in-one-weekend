#include "color.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>

Color color_add(Color c1, Color c2) {
  return (Color){c1.r + c2.r, c1.g + c2.g, c1.b + c2.b};
}

Color color_mul(Color c1, Color c2) {
  return (Color){c1.r * c2.r, c1.g * c2.g, c1.b * c2.b};
}

Color color_lerp(Color c1, Color c2, double t) {
  return (Color){
      (1.0 - t) * c1.r + t * c2.r,
      (1.0 - t) * c1.g + t * c2.g,
      (1.0 - t) * c1.b + t * c2.b,
  };
}

void color_write(FILE *out, Color c, int samples_per_pixel) {
  double scale = 1.0 / samples_per_pixel;
  double r = sqrt(c.r * scale);
  double g = sqrt(c.g * scale);
  double b = sqrt(c.b * scale);

  int ir = (int)(256 * clamp(r, 0.0, 0.999));
  int ig = (int)(256 * clamp(g, 0.0, 0.999));
  int ib = (int)(256 * clamp(b, 0.0, 0.999));

  fprintf(out, "%d %d %d\n", ir, ig, ib);
}

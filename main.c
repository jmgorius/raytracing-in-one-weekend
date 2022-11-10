#include <stdio.h>

#include "color.h"

int main(void) {
  const int image_width = 256;
  const int image_height = 256;

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      color pixel_color = {
          .r = (double)(i) / (image_width - 1),
          .g = (double)(j) / (image_height - 1),
          .b = 0.25,
      };
      color_write(stdout, pixel_color);
    }
  }

  fprintf(stderr, "\nDone.\n");

  return 0;
}

#include "color.c"
#include "point3.c"
#include "vec3.c"

#include "color.h"
#include <stdio.h>

color color_lerp(color c1, color c2, double t) {
  return (color){
      (1.0 - t) * c1.r + t * c2.r,
      (1.0 - t) * c1.g + t * c2.g,
      (1.0 - t) * c1.b + t * c2.b,
  };
}

void color_write(FILE *out, color c) {
  fprintf(out, "%d %d %d\n", (int)(255.999 * c.r), (int)(255.999 * c.g),
          (int)(255.999 * c.b));
}

#include "color.h"
#include <stdio.h>

Color color_lerp(Color c1, Color c2, double t) {
  return (Color){
      (1.0 - t) * c1.r + t * c2.r,
      (1.0 - t) * c1.g + t * c2.g,
      (1.0 - t) * c1.b + t * c2.b,
  };
}

void color_write(FILE *out, Color c) {
  fprintf(out, "%d %d %d\n", (int)(255.999 * c.r), (int)(255.999 * c.g),
          (int)(255.999 * c.b));
}

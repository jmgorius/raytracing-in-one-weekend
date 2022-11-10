#include "color.h"
#include <stdio.h>

void color_write(FILE *out, color c) {
  fprintf(out, "%d %d %d\n", (int)(255.999 * c.r), (int)(255.999 * c.g),
          (int)(255.999 * c.b));
}

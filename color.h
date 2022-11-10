#ifndef INCLUDED_COLOR_H
#define INCLUDED_COLOR_H

#include <stdio.h>

typedef struct Color {
  double r, g, b;
} Color;

Color color_add(Color c1, Color c2);

Color color_lerp(Color c1, Color c2, double t);

void color_write(FILE *out, Color c, int samples_per_pixel);

#endif /* INCLUDED_COLOR_H */

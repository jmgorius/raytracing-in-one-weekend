#ifndef INCLUDED_COLOR_H
#define INCLUDED_COLOR_H

#include <stdio.h>

typedef struct color {
  double r, g, b;
} color;

color color_lerp(color c1, color c2, double t);

void color_write(FILE *out, color c);

#endif /* INCLUDED_COLOR_H */

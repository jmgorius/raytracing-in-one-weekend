#ifndef INCLUDED_COLOR_H
#define INCLUDED_COLOR_H

#include <stdio.h>

typedef struct color {
  double r, g, b;
} color;

void color_write(FILE *out, color c);

#endif /* INCLUDED_COLOR_H */

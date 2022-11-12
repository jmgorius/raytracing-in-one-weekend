#include "utils.h"

#include <math.h>
#include <stdlib.h>

double degrees_to_radians(double degrees) { return degrees * M_PI / 180.0; }

double random_double(void) { return rand() / (RAND_MAX + 1.0); }

double random_double_in_range(double min, double max) {
  return min + (max - min) * random_double();
}

int random_int_in_range(int min, int max) {
  return (int)(random_double_in_range(min, max + 1));
}

double clamp(double x, double min, double max) {
  if (x < min)
    return min;
  if (x > max)
    return max;
  return x;
}

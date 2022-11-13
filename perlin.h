#ifndef INCLUDED_PERLIN_H
#define INCLUDED_PERLIN_H

#include "arena.h"
#include "point3.h"

#define PERLIN_DEFAULT_POINT_COUNT 256
#define PERLIN_DEFAULT_TURBULENCE_DEPTH 7

typedef struct PerlinData PerlinData;

PerlinData *perlin_init(int point_count, Arena *arena);

double perlin_noise(const PerlinData *data, Point3 p);
double perlin_turbulence(const PerlinData *data, Point3 p, int depth);

#endif /* INCLUDED_PERLIN_H */

#ifndef INCLUDED_PERLIN_H
#define INCLUDED_PERLIN_H

#include "arena.h"
#include "point3.h"

#define PERLIN_DEFAULT_POINT_COUNT 256

typedef struct PerlinData PerlinData;

PerlinData *perlin_init(int point_count, Arena *arena);
double perlin_noise(const PerlinData *data, Point3 p);

#endif /* INCLUDED_PERLIN_H */

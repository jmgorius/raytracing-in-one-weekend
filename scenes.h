#ifndef INCLUDED_SCENES_H
#define INCLUDED_SCENES_H

#include "arena.h"
#include "hittable.h"

#define SCENE_RANDOM 0
#define SCENE_TWO_SPHERES 1
#define SCENE_TWO_PERLIN_SPHERES 2
#define SCENE_EARTH 3
#define SCENE_SIMPLE_LIGHT 4
#define SCENE_CORNELL_BOX 5
#define SCENE_CORNELL_BOX_SMOKE 6
#define SCENE_COMPLEX 7

const Hittable *random_scene(Arena *arena);
const Hittable *two_spheres(Arena *arena);
const Hittable *two_perlin_spheres(Arena *arena);
const Hittable *earth(Arena *arena);
const Hittable *simple_light(Arena *arena);
const Hittable *cornell_box(Arena *arena);
const Hittable *cornell_box_smoke(Arena *arena);
const Hittable *complex_scene(Arena *arena);

#endif /* INCLUDED_SCENES_H */

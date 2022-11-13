#include "texture.h"
#include "arena.h"
#include "color.h"
#include "perlin.h"

#include <math.h>

Texture *texture_create_solid_color(Color value, Arena *arena) {
  Texture *result = arena_alloc(arena, sizeof(Texture));
  result->type = TEXTURE_SOLID_COLOR;
  result->solid_color.value = value;
  return result;
}

Texture *texture_create_checker(const Texture *odd, const Texture *even,
                                double size, Arena *arena) {
  Texture *result = arena_alloc(arena, sizeof(Texture));
  result->type = TEXTURE_CHECKER;
  result->checker.odd = odd;
  result->checker.even = even;
  result->checker.size = size;
  return result;
}

Texture *texture_create_perlin_noise(double scale, int point_count,
                                     Arena *arena) {
  Texture *result = arena_alloc(arena, sizeof(Texture));
  result->type = TEXTURE_PERLIN_NOISE;
  result->noise.data = perlin_init(point_count, arena);
  result->noise.scale = scale;
  return result;
}

Texture *texture_create_checker_solid_color(Color odd, Color even, double size,
                                            Arena *arena) {
  Texture *odd_texture = texture_create_solid_color(odd, arena);
  Texture *even_texture = texture_create_solid_color(even, arena);
  return texture_create_checker(odd_texture, even_texture, size, arena);
}

Color texture_value(const Texture *texture, double u, double v, Point3 p) {
  switch (texture->type) {
  case TEXTURE_SOLID_COLOR:
    return solid_color_value(&texture->solid_color);
  case TEXTURE_CHECKER:
    return checker_value(&texture->checker, u, v, p);
  case TEXTURE_PERLIN_NOISE:
    return perlin_value(&texture->noise, p);
  }
  return (Color){0.0, 0.0, 0.0};
}

Color solid_color_value(const SolidColor *solid_color) {
  return solid_color->value;
}

Color checker_value(const CheckerTexture *checker, double u, double v,
                    Point3 p) {
  double sines = sin(checker->size * p.x) * sin(checker->size * p.y) *
                 sin(checker->size * p.z);
  if (sines < 0.0)
    return texture_value(checker->odd, u, v, p);
  return texture_value(checker->even, u, v, p);
}

Color perlin_value(const PerlinNoiseTexture *perlin, Point3 p) {
  double noise_value =
      0.5 * (1.0 + perlin_noise(perlin->data, (Point3){perlin->scale * p.x,
                                                       perlin->scale * p.y,
                                                       perlin->scale * p.z}));
  return color_mul_const(noise_value, (Color){1.0, 1.0, 1.0});
}

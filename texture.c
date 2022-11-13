#include "texture.h"
#include "arena.h"

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

Texture *texture_create_checker_solid_color(Color odd, Color even, double size,
                                            Arena *arena) {
  Texture *odd_texture = texture_create_solid_color(odd, arena);
  Texture *even_texture = texture_create_solid_color(even, arena);
  return texture_create_checker(odd_texture, even_texture, size, arena);
}

Color texture_value(const Texture *texture, double u, double v, Point3 p) {
  switch (texture->type) {
  case TEXTURE_SOLID_COLOR:
    return solid_color_value(&texture->solid_color, u, v, p);
  case TEXTURE_CHECKER:
    return checker_value(&texture->checker, u, v, p);
  }
  return (Color){0.0, 0.0, 0.0};
}

Color solid_color_value(const SolidColor *solid_color, double u, double v,
                        Point3 p) {
  (void)u, (void)v, (void)p;
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

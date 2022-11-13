#include "texture.h"
#include "arena.h"
#include "color.h"
#include "external/stb_image.h"
#include "perlin.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>

#define DEBUG_COLOR                                                            \
  (Color) { 1.0, 0.0, 1.0 }

static const int image_bytes_per_pixel = 3;

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

Texture *texture_create_perlin_noise(double scale, int point_count,
                                     Arena *arena) {
  Texture *result = arena_alloc(arena, sizeof(Texture));
  result->type = TEXTURE_PERLIN_NOISE;
  result->noise.data = perlin_init(point_count, arena);
  result->noise.scale = scale;
  return result;
}

Texture *texture_create_image(const char *filename, Arena *arena) {
  Texture *result = arena_alloc(arena, sizeof(Texture));
  result->type = TEXTURE_IMAGE;

  int components_per_pixel = image_bytes_per_pixel;
  result->image.data =
      stbi_load(filename, &result->image.width, &result->image.height,
                &components_per_pixel, components_per_pixel);
  if (!result->image.data) {
    fprintf(stderr, "Failed to load texture from image file %s\n", filename);
    result->image.width = result->image.height = 0;
  }

  result->image.bytes_per_scanline =
      image_bytes_per_pixel * result->image.width;

  return result;
}

static Color image_value(const ImageTexture *image, double u, double v) {
  if (!image->data)
    return DEBUG_COLOR;

  u = clamp(u, 0.0, 1.0);
  v = 1.0 - clamp(v, 0.0, 1.0);

  int i = u * image->width;
  int j = v * image->height;

  if (i >= image->width)
    i = image->width - 1;
  if (j >= image->height)
    j = image->height - 1;

  const double color_scale = 1.0 / 255.0;
  const unsigned char *pixel =
      image->data + j * image->bytes_per_scanline + i * image_bytes_per_pixel;

  return (Color){color_scale * pixel[0], color_scale * pixel[1],
                 color_scale * pixel[2]};
}

static Color solid_color_value(const SolidColor *solid_color) {
  return solid_color->value;
}

static Color checker_value(const CheckerTexture *checker, double u, double v,
                           Point3 p) {
  double sines = sin(checker->size * p.x) * sin(checker->size * p.y) *
                 sin(checker->size * p.z);
  if (sines < 0.0)
    return texture_value(checker->odd, u, v, p);
  return texture_value(checker->even, u, v, p);
}

static Color perlin_value(const PerlinNoiseTexture *perlin, Point3 p) {
  double noise_value =
      0.5 *
      (1.0 + sin(perlin->scale * p.z +
                 10 * perlin_turbulence(perlin->data, p,
                                        PERLIN_DEFAULT_TURBULENCE_DEPTH)));
  return color_mul_const(noise_value, (Color){1.0, 1.0, 1.0});
}

Color texture_value(const Texture *texture, double u, double v, Point3 p) {
  switch (texture->type) {
  case TEXTURE_SOLID_COLOR:
    return solid_color_value(&texture->solid_color);
  case TEXTURE_CHECKER:
    return checker_value(&texture->checker, u, v, p);
  case TEXTURE_PERLIN_NOISE:
    return perlin_value(&texture->noise, p);
  case TEXTURE_IMAGE:
    return image_value(&texture->image, u, v);
  }
  return DEBUG_COLOR;
}

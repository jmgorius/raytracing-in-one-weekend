#ifndef INCLUDED_TEXTURE_H
#define INCLUDED_TEXTURE_H

#include "arena.h"
#include "color.h"
#include "perlin.h"
#include "point3.h"

typedef enum TextureType {
  TEXTURE_SOLID_COLOR,
  TEXTURE_CHECKER,
  TEXTURE_PERLIN_NOISE,
  TEXTURE_IMAGE,
} TextureType;

typedef struct Texture Texture;

typedef struct SolidColor {
  Color value;
} SolidColor;

typedef struct CheckerTexture {
  const Texture *odd;
  const Texture *even;
  double size;
} CheckerTexture;

typedef struct PerlinNoiseTexture {
  const PerlinData *data;
  double scale;
} PerlinNoiseTexture;

typedef struct ImageTexture {
  unsigned char *data;
  int width, height;
  int bytes_per_scanline;
} ImageTexture;

struct Texture {
  TextureType type;
  union {
    SolidColor solid_color;
    CheckerTexture checker;
    PerlinNoiseTexture noise;
    ImageTexture image;
  };
};

Texture *texture_create_solid_color(Color value, Arena *arena);
Texture *texture_create_checker(const Texture *odd, const Texture *even,
                                double size, Arena *arena);
Texture *texture_create_checker_solid_color(Color odd, Color even, double size,
                                            Arena *arena);
Texture *texture_create_perlin_noise(double scale, int point_count,
                                     Arena *arena);
Texture *texture_create_image(const char *filename, Arena *arena);

Color texture_value(const Texture *texture, double u, double v, Point3 p);

#endif /* INCLUDED_TEXTURE_H */

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

struct Texture {
  TextureType type;
  union {
    SolidColor solid_color;
    CheckerTexture checker;
    PerlinNoiseTexture noise;
  };
};

Texture *texture_create_solid_color(Color value, Arena *arena);
Texture *texture_create_checker(const Texture *odd, const Texture *even,
                                double size, Arena *arena);
Texture *texture_create_checker_solid_color(Color odd, Color even, double size,
                                            Arena *arena);
Texture *texture_create_perlin_noise(double scale, int point_count,
                                     Arena *arena);

Color texture_value(const Texture *texture, double u, double v, Point3 p);

Color solid_color_value(const SolidColor *solid_color);

Color checker_value(const CheckerTexture *checker, double u, double v,
                    Point3 p);

Color perlin_value(const PerlinNoiseTexture *perlin, Point3 p);

#endif /* INCLUDED_TEXTURE_H */

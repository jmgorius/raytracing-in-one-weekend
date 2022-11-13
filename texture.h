#ifndef INCLUDED_TEXTURE_H
#define INCLUDED_TEXTURE_H

#include "arena.h"
#include "color.h"
#include "point3.h"

typedef enum TextureType {
  TEXTURE_SOLID_COLOR,
  TEXTURE_CHECKER,
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

struct Texture {
  TextureType type;
  union {
    SolidColor solid_color;
    CheckerTexture checker;
  };
};

Texture *texture_create_solid_color(Color value, Arena *arena);
Texture *texture_create_checker(const Texture *odd, const Texture *even,
                                double size, Arena *arena);
Texture *texture_create_checker_solid_color(Color odd, Color even, double size,
                                            Arena *arena);

Color texture_value(const Texture *texture, double u, double v, Point3 p);

Color solid_color_value(const SolidColor *solid_color, double u, double v,
                        Point3 p);

Color checker_value(const CheckerTexture *checker, double u, double v,
                    Point3 p);

#endif /* INCLUDED_TEXTURE_H */

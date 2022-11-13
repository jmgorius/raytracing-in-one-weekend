#ifndef INCLUDED_MATERIAL_H
#define INCLUDED_MATERIAL_H

#include "arena.h"
#include "ray.h"
#include "texture.h"

#include <stdbool.h>

struct HitRecord;

typedef enum MaterialType {
  MATERIAL_LAMBERTIAN,
  MATERIAL_METAL,
  MATERIAL_DIELECTRIC,
} MaterialType;

typedef struct Lambertian {
  const Texture *albedo;
} Lambertian;

typedef struct Metal {
  const Texture *albedo;
  double fuzziness;
} Metal;

typedef struct Dielectric {
  double eta;
} Dielectric;

typedef struct Material {
  MaterialType type;
  union {
    Lambertian lambertian;
    Metal metal;
    Dielectric dielectric;
  };
} Material;

bool material_scatter(const Material *material, Ray r,
                      const struct HitRecord *record, Color *attenuation,
                      Ray *scattered);

Material *material_create_lambertian(const Texture *albedo, Arena *arena);
Material *material_create_lambertian_color(Color albedo, Arena *arena);
Material *material_create_metal(const Texture *albedo, double fuzziness,
                                Arena *arena);
Material *material_create_metal_color(Color albedo, double fuzziness,
                                      Arena *arena);
Material *material_create_dielectric(double eta, Arena *arena);

#endif /* INCLUDED_MATERIAL_H */

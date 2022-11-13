#ifndef INCLUDED_MATERIAL_H
#define INCLUDED_MATERIAL_H

#include "arena.h"
#include "color.h"
#include "ray.h"

#include <stdbool.h>

struct HitRecord;

typedef enum MaterialType {
  MATERIAL_LAMBERTIAN,
  MATERIAL_METAL,
  MATERIAL_DIELECTRIC,
} MaterialType;

typedef struct Lambertian {
  Color albedo;
} Lambertian;

typedef struct Metal {
  Color albedo;
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

Material *material_create_lambertian(Color albedo, Arena *arena);
bool lambertian_scatter(const Lambertian *lambertian, Ray r,
                        const struct HitRecord *record, Color *attenuation,
                        Ray *scattered);

Material *material_create_metal(Color albedo, double fuzziness, Arena *arena);
bool metal_scatter(const Metal *metal, Ray r, const struct HitRecord *record,
                   Color *attenuation, Ray *scattered);

Material *material_create_dielectric(double eta, Arena *arena);
bool dielectric_scatter(const Dielectric *dielectric, Ray r,
                        const struct HitRecord *record, Color *attenuation,
                        Ray *scattered);

#endif /* INCLUDED_MATERIAL_H */

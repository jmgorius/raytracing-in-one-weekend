#ifndef INCLUDED_MATERIAL_H
#define INCLUDED_MATERIAL_H

#include "color.h"
#include "ray.h"

#include <stdbool.h>

struct HitRecord;

typedef enum MaterialType {
  MATERIAL_LAMBERTIAN,
  MATERIAL_METAL,
  MATERIAL_DIELECTRIC,
} MaterialType;

typedef struct Material {
  MaterialType type;
} Material;

bool material_scatter(const Material *material, Ray r,
                      const struct HitRecord *record, Color *attenuation,
                      Ray *scattered);

typedef struct Lambertian {
  MaterialType type;
  Color albedo;
} Lambertian;

bool lambertian_scatter(const Lambertian *lambertian, Ray r,
                        const struct HitRecord *record, Color *attenuation,
                        Ray *scattered);

typedef struct Metal {
  MaterialType type;
  Color albedo;
  double fuzziness;
} Metal;

bool metal_scatter(const Metal *metal, Ray r, const struct HitRecord *record,
                   Color *attenuation, Ray *scattered);

typedef struct Dielectric {
  MaterialType type;
  double eta;
} Dielectric;

bool dielectric_scatter(const Dielectric *dielectric, Ray r,
                        const struct HitRecord *record, Color *attenuation,
                        Ray *scattered);

#endif /* INCLUDED_MATERIAL_H */

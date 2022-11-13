#include "material.h"
#include "arena.h"
#include "hittable.h"
#include "utils.h"
#include "vec3.h"

#include <assert.h>
#include <math.h>

static bool lambertian_scatter(const Lambertian *lambertian, Ray r,
                               const HitRecord *record, Color *attenuation,
                               Ray *scattered) {
  (void)r;
  Vec3 scatter_direction = vec3_add(record->normal, vec3_random_unit_vector());

  /* Catch degenerate scatter direction */
  if (vec3_is_near_zero(scatter_direction))
    scatter_direction = record->normal;

  *scattered = (Ray){record->p, scatter_direction, r.time};
  *attenuation =
      texture_value(lambertian->albedo, record->u, record->v, record->p);
  return true;
}

static bool metal_scatter(const Metal *metal, Ray r,
                          const struct HitRecord *record, Color *attenuation,
                          Ray *scattered) {
  Vec3 reflected = vec3_reflect(vec3_normalize(r.direction), record->normal);
  assert(metal->fuzziness <= 1);
  *scattered = (Ray){
      record->p,
      vec3_add(reflected,
               vec3_mul(metal->fuzziness, vec3_random_in_unit_sphere())),
      r.time,
  };
  *attenuation = texture_value(metal->albedo, record->u, record->v, record->p);
  return vec3_dot(scattered->direction, record->normal) > 0;
}

static double schlick_reflectance(double cosine, double eta) {
  double r0 = (1 - eta) / (1 + eta);
  r0 *= r0;
  return r0 + (1 - r0) * pow(1 - cosine, 5);
}

static bool dielectric_scatter(const Dielectric *dielectric, Ray r,
                               const struct HitRecord *record,
                               Color *attenuation, Ray *scattered) {
  *attenuation = (Color){1.0, 1.0, 1.0};
  double refraction_ratio =
      record->front_face ? (1.0 / dielectric->eta) : dielectric->eta;

  Vec3 unit_direction = vec3_normalize(r.direction);
  double cos_theta =
      fmin(vec3_dot(vec3_neg(unit_direction), record->normal), 1.0);
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  bool cannot_refract = refraction_ratio * sin_theta > 1.0;

  Vec3 direction;
  if (cannot_refract ||
      schlick_reflectance(cos_theta, refraction_ratio) > random_double())
    direction = vec3_reflect(unit_direction, record->normal);
  else
    direction = vec3_refract(unit_direction, record->normal, refraction_ratio);

  *scattered = (Ray){record->p, direction, r.time};
  return true;
}

bool material_scatter(const Material *material, Ray r,
                      const struct HitRecord *record, Color *attenuation,
                      Ray *scattered) {
  switch (material->type) {
  case MATERIAL_LAMBERTIAN:
    return lambertian_scatter(&material->lambertian, r, record, attenuation,
                              scattered);
  case MATERIAL_METAL:
    return metal_scatter(&material->metal, r, record, attenuation, scattered);
  case MATERIAL_DIELECTRIC:
    return dielectric_scatter(&material->dielectric, r, record, attenuation,
                              scattered);
  }
  return false;
}

Material *material_create_lambertian(const Texture *albedo, Arena *arena) {
  Material *result = arena_alloc(arena, sizeof(Material));
  result->type = MATERIAL_LAMBERTIAN;
  result->lambertian.albedo = albedo;
  return result;
}

Material *material_create_lambertian_color(Color albedo, Arena *arena) {
  return material_create_lambertian(texture_create_solid_color(albedo, arena),
                                    arena);
}

Material *material_create_metal(const Texture *albedo, double fuzziness,
                                Arena *arena) {
  Material *result = arena_alloc(arena, sizeof(Material));
  result->type = MATERIAL_METAL;
  result->metal.albedo = albedo;
  result->metal.fuzziness = fuzziness;
  return result;
}

Material *material_create_metal_color(Color albedo, double fuzziness,
                                      Arena *arena) {
  return material_create_metal(texture_create_solid_color(albedo, arena),
                               fuzziness, arena);
}

Material *material_create_dielectric(double eta, Arena *arena) {
  Material *result = arena_alloc(arena, sizeof(Material));
  result->type = MATERIAL_DIELECTRIC;
  result->dielectric.eta = eta;
  return result;
}

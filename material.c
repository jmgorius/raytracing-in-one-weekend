#include "material.h"
#include "hittable.h"
#include "vec3.h"

#include <assert.h>

bool material_scatter(const Material *material, Ray r,
                      const struct HitRecord *record, Color *attenuation,
                      Ray *scattered) {
  switch (material->type) {
  case MATERIAL_LAMBERTIAN:
    return lambertian_scatter((const Lambertian *)material, r, record,
                              attenuation, scattered);
  case MATERIAL_METAL:
    return metal_scatter((const Metal *)material, r, record, attenuation,
                         scattered);
  }
  return false;
}

bool lambertian_scatter(const Lambertian *lambertian, Ray r,
                        const HitRecord *record, Color *attenuation,
                        Ray *scattered) {
  (void)r;
  Vec3 scatter_direction = vec3_add(record->normal, vec3_random_unit_vector());

  /* Catch degenerate scatter direction */
  if (vec3_is_near_zero(scatter_direction))
    scatter_direction = record->normal;

  *scattered = (Ray){record->p, scatter_direction};
  *attenuation = lambertian->albedo;
  return true;
}

bool metal_scatter(const Metal *metal, Ray r, const struct HitRecord *record,
                   Color *attenuation, Ray *scattered) {
  Vec3 reflected = vec3_reflect(vec3_normalize(r.direction), record->normal);
  assert(metal->fuzziness <= 1);
  *scattered = (Ray){
      record->p,
      vec3_add(reflected,
               vec3_mul(metal->fuzziness, vec3_random_in_unit_sphere())),
  };
  *attenuation = metal->albedo;
  return vec3_dot(scattered->direction, record->normal) > 0;
}

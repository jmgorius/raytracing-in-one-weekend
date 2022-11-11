#include "camera.h"
#include "point3.h"
#include "utils.h"
#include "vec3.h"

#include <math.h>

void camera_init(Camera *camera, Point3 look_from, Point3 look_at, Vec3 up,
                 double vertical_fov, double aspect_ratio) {
  double theta = degrees_to_radians(vertical_fov);
  double h = tan(theta / 2.0);
  double viewport_height = 2.0 * h;
  double viewport_width = aspect_ratio * viewport_height;

  Vec3 w = vec3_normalize(point3_sub(look_from, look_at));
  Vec3 u = vec3_normalize(vec3_cross(up, w));
  Vec3 v = vec3_cross(w, u);

  camera->origin = look_from;
  camera->horizontal = vec3_mul(viewport_width, u);
  camera->vertical = vec3_mul(viewport_height, v);
  Vec3 offset =
      vec3_add(vec3_div(camera->horizontal, 2), vec3_div(camera->vertical, 2));
  offset = vec3_add(offset, w);
  camera->lower_left_corner = point3_add(camera->origin, vec3_neg(offset));
}

Ray camera_get_ray(const Camera *camera, double s, double t) {
  Point3 screen_point = point3_add(
      camera->lower_left_corner,
      vec3_add(vec3_mul(s, camera->horizontal), vec3_mul(t, camera->vertical)));
  Vec3 direction = point3_sub(screen_point, camera->origin);
  return (Ray){camera->origin, direction};
}

#include "camera.h"

void camera_init(Camera *camera, double aspect_ratio) {
  double viewport_height = 2;
  double viewport_width = aspect_ratio * viewport_height;
  double focal_length = 1.0;

  camera->origin = (Point3){0, 0, 0};
  camera->horizontal = (Vec3){viewport_width, 0, 0};
  camera->vertical = (Vec3){0, viewport_height, 0};
  Vec3 offset =
      vec3_add(vec3_div(camera->horizontal, 2), vec3_div(camera->vertical, 2));
  offset = vec3_add(offset, (Vec3){0, 0, focal_length});
  camera->lower_left_corner = point3_add(camera->origin, vec3_neg(offset));
}

Ray camera_get_ray(const Camera *camera, double u, double v) {
  Point3 screen_point = point3_add(
      camera->lower_left_corner,
      vec3_add(vec3_mul(u, camera->horizontal), vec3_mul(v, camera->vertical)));
  Vec3 direction = point3_sub(screen_point, camera->origin);
  return (Ray){camera->origin, direction};
}

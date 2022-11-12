#include "camera.h"
#include "point3.h"
#include "utils.h"
#include "vec3.h"

#include <math.h>

void camera_init(Camera *camera, Point3 look_from, Point3 look_at, Vec3 up,
                 double vertical_fov, double aspect_ratio, double aperture,
                 double focus_distance, double shutter_open,
                 double shutter_close) {
  double theta = degrees_to_radians(vertical_fov);
  double h = tan(theta / 2.0);
  double viewport_height = 2.0 * h;
  double viewport_width = aspect_ratio * viewport_height;

  camera->w = vec3_normalize(point3_sub(look_from, look_at));
  camera->u = vec3_normalize(vec3_cross(up, camera->w));
  camera->v = vec3_cross(camera->w, camera->u);

  camera->origin = look_from;
  camera->horizontal = vec3_mul(focus_distance * viewport_width, camera->u);
  camera->vertical = vec3_mul(focus_distance * viewport_height, camera->v);
  Vec3 offset =
      vec3_add(vec3_div(camera->horizontal, 2), vec3_div(camera->vertical, 2));
  offset = vec3_add(offset, vec3_mul(focus_distance, camera->w));
  camera->lower_left_corner = point3_add(camera->origin, vec3_neg(offset));

  camera->lens_radius = aperture / 2.0;
  camera->shutter_open = shutter_open;
  camera->shutter_close = shutter_close;
}

Ray camera_get_ray(const Camera *camera, double s, double t) {
  Vec3 rd = vec3_mul(camera->lens_radius, vec3_random_in_unit_disk());
  Vec3 offset = vec3_add(vec3_mul(rd.x, camera->u), vec3_mul(rd.y, camera->v));

  Point3 screen_point = point3_add(
      camera->lower_left_corner,
      vec3_add(vec3_mul(s, camera->horizontal), vec3_mul(t, camera->vertical)));
  Vec3 direction = point3_sub(screen_point, camera->origin);
  return (Ray){
      .origin = point3_add(camera->origin, offset),
      .direction = vec3_sub(direction, offset),
      .time =
          random_double_in_range(camera->shutter_open, camera->shutter_close),
  };
}

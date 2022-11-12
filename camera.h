#ifndef INCLUDED_CAMERA_H
#define INCLUDED_CAMERA_H

#include "point3.h"
#include "ray.h"
#include "vec3.h"

typedef struct Camera {
  Point3 origin;
  Point3 lower_left_corner;
  Vec3 horizontal, vertical;
  Vec3 u, v, w;
  double lens_radius;
  double shutter_open, shutter_close;
} Camera;

void camera_init(Camera *camera, Point3 look_from, Point3 look_at, Vec3 up,
                 double vertical_fov, double aspect_ratio, double aperture,
                 double focus_distance, double shutter_open,
                 double shutter_close);

Ray camera_get_ray(const Camera *camera, double s, double t);

#endif /* INCLUDED_CAMERA_H */

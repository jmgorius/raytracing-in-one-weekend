#ifndef INCLUDED_CAMERA_H
#define INCLUDED_CAMERA_H

#include "point3.h"
#include "ray.h"
#include "vec3.h"

typedef struct Camera {
  Point3 origin;
  Point3 lower_left_corner;
  Vec3 horizontal;
  Vec3 vertical;
} Camera;

void camera_init(Camera *camera, double aspect_ratio);

Ray camera_get_ray(const Camera *camera, double u, double v);

#endif /* INCLUDED_CAMERA_H */

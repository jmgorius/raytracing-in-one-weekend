#include "point3.h"

Point3 point3_add(Point3 p, Vec3 v) {
  return (Point3){p.x + v.x, p.y + v.y, p.z + v.z};
}

Vec3 point3_sub(Point3 p1, Point3 p2) {
  return (Vec3){p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}

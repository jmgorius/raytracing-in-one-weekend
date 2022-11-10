#include "point3.h"

point3 point3_add(point3 p, vec3 v) {
  return (point3){p.x + v.x, p.y + v.y, p.z + v.z};
}

vec3 point3_sub(point3 p1, point3 p2) {
  return (vec3){p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}

#include "perlin.h"
#include "arena.h"
#include "utils.h"
#include "vec3.h"

#include <math.h>

struct PerlinData {
  Vec3 *random_vectors;
  int *perm_x;
  int *perm_y;
  int *perm_z;
};

static void permute(int *p, int n) {
  for (int i = n - 1; i > 0; --i) {
    int target = random_int_in_range(0, i);
    int tmp = p[i];
    p[i] = p[target];
    p[target] = tmp;
  }
}

static int *perlin_generate_perm(int point_count, Arena *arena) {
  int *p = arena_alloc(arena, point_count * sizeof(int));

  for (int i = 0; i < point_count; ++i)
    p[i] = i;
  permute(p, point_count);
  return p;
}

PerlinData *perlin_init(int point_count, Arena *arena) {
  PerlinData *data = arena_alloc(arena, sizeof(PerlinData));
  data->random_vectors = arena_alloc(arena, point_count * sizeof(Vec3));
  for (int i = 0; i < point_count; ++i)
    data->random_vectors[i] = vec3_normalize(vec3_random_in_range(-1.0, 1.0));

  data->perm_x = perlin_generate_perm(point_count, arena);
  data->perm_y = perlin_generate_perm(point_count, arena);
  data->perm_z = perlin_generate_perm(point_count, arena);

  return data;
}

static double perlin_interp(const Vec3 c[2][2][2], double u, double v,
                            double w) {
  double uu = u * u * (3 - 2 * u);
  double vv = v * v * (3 - 2 * v);
  double ww = w * w * (3 - 2 * w);
  double accum = 0.0;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        Vec3 weight_v = {u - i, v - j, w - k};
        accum += (i * uu + (1 - i) * (1 - uu)) * (j * vv + (1 - j) * (1 - vv)) *
                 (k * ww + (1 - k) * (1 - ww)) * vec3_dot(c[i][j][k], weight_v);
      }
    }
  }
  return accum;
}

double perlin_noise(const PerlinData *data, Point3 p) {
  double u = p.x - floor(p.x);
  double v = p.y - floor(p.y);
  double w = p.z - floor(p.z);

  int i = (int)floor(p.x);
  int j = (int)floor(p.y);
  int k = (int)floor(p.z);

  Vec3 c[2][2][2] = {0};
  for (int di = 0; di < 2; ++di) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int dk = 0; dk < 2; ++dk) {
        c[di][dj][dk] = data->random_vectors[data->perm_x[(i + di) & 0xff] ^
                                             data->perm_y[(j + dj) & 0xff] ^
                                             data->perm_z[(k + dk) & 0xff]];
      }
    }
  }

  return perlin_interp(c, u, v, w);
}

#include "scenes.h"
#include "utils.h"

const Hittable *random_scene(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  const Texture *checker = texture_create_checker_solid_color(
      (Color){0.2, 0.3, 0.1}, (Color){0.9, 0.9, 0.9}, 10.0, arena);
  const Material *ground_material = material_create_lambertian(checker, arena);
  const Hittable *ground_sphere = hittable_create_sphere(
      (Point3){0.0, -1000.0, 0.0}, 1000.0, ground_material, arena);
  hittable_list_add(&world->list, ground_sphere, arena);

  const Material *glass = material_create_dielectric(1.5, arena);

  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      double choose_material = random_double();
      Point3 center = {a + 0.9 * random_double(), 0.2,
                       b + 0.9 * random_double()};

      if (vec3_length(point3_sub(center, (Point3){4.0, 0.2, 0.0})) > 0.9) {
        if (choose_material < 0.8) {
          const Material *material = material_create_lambertian_color(
              color_mul(color_random(), color_random()), arena);
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, material, arena);
          hittable_list_add(&world->list, sphere, arena);
        } else if (choose_material < 0.95) {
          const Material *material = material_create_metal_color(
              color_random_in_range(0.5, 1), random_double_in_range(0.5, 1.0),
              arena);
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, material, arena);
          hittable_list_add(&world->list, sphere, arena);
        } else {
          const Hittable *sphere =
              hittable_create_sphere(center, 0.2, glass, arena);
          hittable_list_add(&world->list, sphere, arena);
        }
      }
    }
  }

  const Material *lambertian =
      material_create_lambertian_color((Color){0.4, 0.2, 0.1}, arena);
  const Material *metal =
      material_create_metal_color((Color){0.7, 0.6, 0.5}, 0.0, arena);

  const Hittable *sphere1 =
      hittable_create_sphere((Point3){0.0, 1.0, 0.0}, 1.0, glass, arena);
  hittable_list_add(&world->list, sphere1, arena);
  const Hittable *sphere2 =
      hittable_create_sphere((Point3){-4.0, 1.0, 0.0}, 1.0, lambertian, arena);
  hittable_list_add(&world->list, sphere2, arena);
  const Hittable *sphere3 =
      hittable_create_sphere((Point3){4.0, 1.0, 0.0}, 1.0, metal, arena);
  hittable_list_add(&world->list, sphere3, arena);

  Hittable *bvh_root = hittable_create_bvh_node(
      world->list.objects, 0, world->list.size, 0.0, 1.0, arena);
  return bvh_root;
}

const Hittable *two_spheres(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *checker = texture_create_checker_solid_color(
      (Color){0.2, 0.3, 0.1}, (Color){0.9, 0.9, 0.9}, 10.0, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -10.0, 0.0}, 10.0,
                             material_create_lambertian(checker, arena), arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 10.0, 0.0}, 10.0,
                             material_create_lambertian(checker, arena), arena),
      arena);

  return world;
}

const Hittable *two_perlin_spheres(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *perlin_texture =
      texture_create_perlin_noise(4.0, PERLIN_DEFAULT_POINT_COUNT, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -1000.0, 0.0}, 1000.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 2.0, 0.0}, 2.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);

  return world;
}

const Hittable *earth(Arena *arena) {
  Texture *earth_texture = texture_create_image("assets/earthmap.jpg", arena);
  Material *earth_surface = material_create_lambertian(earth_texture, arena);
  Hittable *globe = hittable_create_sphere((Point3){0.0, 0.0, 0.0}, 2.0,
                                           earth_surface, arena);
  return globe;
}

const Hittable *simple_light(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Texture *perlin_texture =
      texture_create_perlin_noise(4.0, PERLIN_DEFAULT_POINT_COUNT, arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, -1000.0, 0.0}, 1000.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_sphere((Point3){0.0, 2.0, 0.0}, 2.0,
                             material_create_lambertian(perlin_texture, arena),
                             arena),
      arena);

  Material *diffuse_light =
      material_create_diffuse_light_color((Color){4.0, 4.0, 4.0}, arena);
  hittable_list_add(&world->list,
                    hittable_create_xy_rectangle(3.0, 5.0, 1.0, 3.0, -2.0,
                                                 diffuse_light, arena),
                    arena);

  return world;
}

const Hittable *cornell_box(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Material *red =
      material_create_lambertian_color((Color){0.65, 0.05, 0.05}, arena);
  Material *white =
      material_create_lambertian_color((Color){0.73, 0.73, 0.73}, arena);
  Material *green =
      material_create_lambertian_color((Color){0.12, 0.45, 0.15}, arena);
  Material *light =
      material_create_diffuse_light_color((Color){15.0, 15.0, 15.0}, arena);

  hittable_list_add(
      &world->list,
      hittable_create_yz_rectangle(0, 555, 0, 555, 555, green, arena), arena);
  hittable_list_add(&world->list,
                    hittable_create_yz_rectangle(0, 555, 0, 555, 0, red, arena),
                    arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(213, 343, 227, 332, 554, light, arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(0, 555, 0, 555, 0, white, arena), arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(0, 555, 0, 555, 555, white, arena), arena);
  hittable_list_add(
      &world->list,
      hittable_create_xy_rectangle(0, 555, 0, 555, 555, white, arena), arena);

  Hittable *box1 = hittable_create_box((Point3){0, 0, 0},
                                       (Point3){165, 330, 165}, white, arena);
  box1 = hittable_create_y_rotation(box1, 15, arena);
  box1 = hittable_create_translation(box1, (Vec3){265, 0, 295}, arena);
  hittable_list_add(&world->list, box1, arena);

  Hittable *box2 = hittable_create_box((Point3){0, 0, 0},
                                       (Point3){165, 165, 165}, white, arena);
  box2 = hittable_create_y_rotation(box2, -18, arena);
  box2 = hittable_create_translation(box2, (Vec3){130, 0, 65}, arena);
  hittable_list_add(&world->list, box2, arena);

  return world;
}

const Hittable *cornell_box_smoke(Arena *arena) {
  Hittable *world = hittable_create_hittable_list(arena);

  Material *red =
      material_create_lambertian_color((Color){0.65, 0.05, 0.05}, arena);
  Material *white =
      material_create_lambertian_color((Color){0.73, 0.73, 0.73}, arena);
  Material *green =
      material_create_lambertian_color((Color){0.12, 0.45, 0.15}, arena);
  Material *light =
      material_create_diffuse_light_color((Color){7.0, 7.0, 7.0}, arena);

  hittable_list_add(
      &world->list,
      hittable_create_yz_rectangle(0, 555, 0, 555, 555, green, arena), arena);
  hittable_list_add(&world->list,
                    hittable_create_yz_rectangle(0, 555, 0, 555, 0, red, arena),
                    arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(113, 443, 127, 432, 554, light, arena),
      arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(0, 555, 0, 555, 0, white, arena), arena);
  hittable_list_add(
      &world->list,
      hittable_create_xz_rectangle(0, 555, 0, 555, 555, white, arena), arena);
  hittable_list_add(
      &world->list,
      hittable_create_xy_rectangle(0, 555, 0, 555, 555, white, arena), arena);

  Hittable *box1 = hittable_create_box((Point3){0, 0, 0},
                                       (Point3){165, 330, 165}, white, arena);
  box1 = hittable_create_y_rotation(box1, 15, arena);
  box1 = hittable_create_translation(box1, (Vec3){265, 0, 295}, arena);
  hittable_list_add(&world->list,
                    hittable_create_constant_medium_color(
                        box1, 0.01, (Color){0.0, 0.0, 0.0}, arena),
                    arena);

  Hittable *box2 = hittable_create_box((Point3){0, 0, 0},
                                       (Point3){165, 165, 165}, white, arena);
  box2 = hittable_create_y_rotation(box2, -18, arena);
  box2 = hittable_create_translation(box2, (Vec3){130, 0, 65}, arena);
  hittable_list_add(&world->list,
                    hittable_create_constant_medium_color(
                        box2, 0.01, (Color){1.0, 1.0, 1.0}, arena),
                    arena);

  return world;
}

const Hittable *complex_scene(Arena *arena) {
  Hittable *boxes1 = hittable_create_hittable_list(arena);

  Material *ground =
      material_create_lambertian_color((Color){0.48, 0.83, 0.53}, arena);

  const int boxes_per_side = 20;
  for (int i = 0; i < boxes_per_side; ++i) {
    for (int j = 0; j < boxes_per_side; ++j) {
      double w = 100.0;
      double x0 = -1000.0 + i * w;
      double y0 = 0.0;
      double z0 = -1000.0 + j * w;
      double x1 = x0 + w;
      double y1 = random_double_in_range(1, 101);
      double z1 = z0 + w;
      hittable_list_add(&boxes1->list,
                        hittable_create_box((Point3){x0, y0, z0},
                                            (Point3){x1, y1, z1}, ground,
                                            arena),
                        arena);
    }
  }

  Hittable *objects = hittable_create_hittable_list(arena);
  hittable_list_add(&objects->list,
                    hittable_create_bvh_node(boxes1->list.objects, 0,
                                             boxes1->list.size, 0, 1, arena),
                    arena);

  Material *light =
      material_create_diffuse_light_color((Color){7, 7, 7}, arena);
  hittable_list_add(
      &objects->list,
      hittable_create_xz_rectangle(123, 423, 147, 412, 554, light, arena),
      arena);

  hittable_list_add(&objects->list,
                    hittable_create_sphere((Point3){400, 400, 200}, 50,
                                           material_create_lambertian_color(
                                               (Color){0.7, 0.3, 0.1}, arena),
                                           arena),
                    arena);
  hittable_list_add(
      &objects->list,
      hittable_create_sphere((Point3){260, 150, 45}, 50,
                             material_create_dielectric(1.5, arena), arena),
      arena);
  hittable_list_add(
      &objects->list,
      hittable_create_sphere(
          (Point3){0, 150, 145}, 50,
          material_create_metal_color((Color){0.8, 0.8, 0.9}, 1.0, arena),
          arena),
      arena);

  Hittable *boundary =
      hittable_create_sphere((Point3){360, 150, 145}, 70,
                             material_create_dielectric(1.5, arena), arena);
  hittable_list_add(&objects->list, boundary, arena);
  hittable_list_add(&objects->list,
                    hittable_create_constant_medium_color(
                        boundary, 0.2, (Color){0.2, 0.4, 0.9}, arena),
                    arena);
  /* Ambient fog */
  boundary = hittable_create_sphere(
      (Point3){0, 0, 0}, 5000, material_create_dielectric(1.5, arena), arena);
  hittable_list_add(&objects->list,
                    hittable_create_constant_medium_color(
                        boundary, 0.0001, (Color){1.0, 1.0, 1.0}, arena),
                    arena);

  hittable_list_add(
      &objects->list,
      hittable_create_sphere(
          (Point3){400, 200, 400}, 100,
          material_create_lambertian(
              texture_create_image("assets/earthmap.jpg", arena), arena),
          arena),
      arena);
  hittable_list_add(
      &objects->list,
      hittable_create_sphere((Point3){220, 280, 300}, 80,
                             material_create_lambertian(
                                 texture_create_perlin_noise(
                                     0.1, PERLIN_DEFAULT_POINT_COUNT, arena),
                                 arena),
                             arena),
      arena);

  Hittable *boxes2 = hittable_create_hittable_list(arena);
  Material *white =
      material_create_lambertian_color((Color){0.73, 0.73, 0.73}, arena);
  int ns = 1000;
  for (int i = 0; i < ns; ++i) {
    hittable_list_add(&boxes2->list,
                      hittable_create_sphere(point3_random_in_range(0, 165), 10,
                                             white, arena),
                      arena);
  }

  hittable_list_add(
      &objects->list,
      hittable_create_translation(
          hittable_create_y_rotation(
              hittable_create_bvh_node(boxes2->list.objects, 0,
                                       boxes2->list.size, 0, 1, arena),
              15, arena),
          (Vec3){-100, 270, 395}, arena),
      arena);

  return objects;
}

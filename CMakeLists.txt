cmake_minimum_required(VERSION 3.10)

project(raytracing-in-one-weekend C)

add_executable(raytracer main.c
               aabb.c
               arena.c
               camera.c
               color.c
               hittable.c
               material.c
               perlin.c
               point3.c
               ray.c
               scenes.c
               texture.c
               utils.c
               vec3.c

               external/stb_image.c)
target_link_libraries(raytracer PRIVATE m)
target_compile_options(raytracer PRIVATE -march=native)

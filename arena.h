#ifndef INCLUDED_ARENA_H
#define INCLUDED_ARENA_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct Arena {
  unsigned char *buffer;
  size_t buffer_length;
  size_t offset;
} Arena;

void arena_init(Arena *arena, void *buffer, size_t buffer_length);

void *arena_alloc(Arena *arena, size_t size);
void *arena_alloc_align(Arena *arena, size_t size, size_t alignment);
void arena_clear(Arena *arena);

#endif /* INCLUDED_ARENA_H */

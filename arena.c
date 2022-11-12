#include "arena.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

static bool is_power_of_two(uintptr_t x) { return (x & (x - 1)) == 0; }

static uintptr_t align_forward(uintptr_t ptr, size_t alignment) {
  assert(is_power_of_two(alignment));

  uintptr_t modulo = ptr & (alignment - 1); /* = ptr % alignment */
  if (modulo != 0)
    ptr += alignment - modulo;
  return ptr;
}

void arena_init(Arena *arena, void *buffer, size_t buffer_length) {
  arena->buffer = buffer;
  arena->buffer_length = buffer_length;
  arena->offset = 0;
}

#ifndef DEFAULT_ALIGNMENT
#define DEFAULT_ALIGNMENT (2 * sizeof(void *))
#endif

void *arena_alloc(Arena *arena, size_t size) {
  return arena_alloc_align(arena, size, DEFAULT_ALIGNMENT);
}

void *arena_alloc_align(Arena *arena, size_t size, size_t alignment) {
  /* Align the current offset to the specified alignment */
  uintptr_t curr_ptr = (uintptr_t)arena->buffer + arena->offset;
  uintptr_t offset = align_forward(curr_ptr, alignment);
  offset -= (uintptr_t)arena->buffer;

  /* Check to see if there is some memory left */
  if (offset <= arena->buffer_length - size) {
    void *ptr = &arena->buffer[offset];
    arena->offset = offset + size;

    /* Zero out newly allocated memory */
    memset(ptr, 0, size);
    assert(ptr != 0);
    return ptr;
  }

  abort();
}

void arena_clear(Arena *arena) { arena->offset = 0; }

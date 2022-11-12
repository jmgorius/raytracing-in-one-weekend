#!/bin/bash

CC=${CC:-"gcc"}
CFLAGS=${CFLAGS:-"-Wall -Wextra -std=gnu18 -fno-strict-aliasing -O3"}
LDFLAGS=${LDFLAGS:-"-lm"}
MAIN_FILE=${MAIN_FILE:-"main.c"}

# Write out compile_flags.txt for clangd
echo "${CFLAGS}" | tr ' ' '\n' > compile_flags.txt

mkdir -p build

${CC} ${CFLAGS} ${MAIN_FILE} -o build/raytracer ${LDFLAGS}

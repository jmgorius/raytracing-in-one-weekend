#include <stdio.h>

int main(void) {
  const int image_width = 256;
  const int image_height = 256;

  printf("P3\n%u %u\n255\n", image_width, image_height);

  for (int j = image_height - 1; j >= 0; --j) {
    fprintf(stderr, "\rScanlines remaining: %d      ", j);
    for (int i = 0; i < image_width; ++i) {
      double r = (double)(i) / (image_width - 1);
      double g = (double)(j) / (image_height - 1);
      double b = 0.25;

      int ir = (int)(255.999 * r);
      int ig = (int)(255.999 * g);
      int ib = (int)(255.999 * b);

      printf("%d %d %d\n", ir, ig, ib);
    }
  }

  fprintf(stderr, "\nDone.\n");

  return 0;
}

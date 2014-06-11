#include <stdio.h>
#include <stdlib.h>

#define WIDTH 200
#define HEIGHT 200

int main() {
  int i, j;
  printf("ply\nelement vertex %d\nelement face %d\nend_header\n", WIDTH*HEIGHT,
	 (WIDTH-1)*(HEIGHT-1)*2);
  for (i = 0; i < HEIGHT; ++i) {
    for (j = 0; j < WIDTH; ++j) { 
      printf("%f %f %f\n", 0.0, (2*i-HEIGHT)*0.001, (2*j-WIDTH)*0.001);
    }
  }
  for (i = 1; i < HEIGHT; ++i) {
    for (j = 1; j < WIDTH; ++j) { 
      printf("3 %d %d %d\n", (i-1)*WIDTH+j-1, i*WIDTH+j-1, (i-1)*WIDTH+j);
      printf("3 %d %d %d\n", (i-1)*WIDTH+j, i*WIDTH+j-1, i*WIDTH+j);
    }
  }
}

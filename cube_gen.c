#include <stdio.h>
#include <stdlib.h>

#define WIDTH 30
#define HEIGHT 30
#define DEPTH 30

int main() {
  int i, j, k, v;
  printf("plt\nelement vertex %d\nelement tet %d\nend_header\n",
         WIDTH*HEIGHT*DEPTH, 5 * (WIDTH-1) * (HEIGHT-1) * (DEPTH-1));
  for (i = 0; i < WIDTH; ++i)
    for (j = 0; j < HEIGHT; ++j)
      for (k = 0; k < DEPTH; ++k) {
        printf("%f %f %f\n", 
               (2*i-WIDTH) * 0.01,
               (2*j-HEIGHT) * 0.01,
               (2*k-DEPTH) * 0.01);
      }

  for (i = 1; i < WIDTH; ++i)
    for (j = 1; j < HEIGHT; ++j)
      for (k = 1; k < DEPTH; ++k) {
        v = i * HEIGHT * DEPTH + j * DEPTH + k;
        printf("4 %d %d %d %d\n", v, v-HEIGHT*DEPTH, v-DEPTH, v-1);
        /*
        printf("4 %d %d %d %d\n", v-(HEIGHT*DEPTH+DEPTH), v-HEIGHT*DEPTH, v-DEPTH, v-(1+HEIGHT*DEPTH+DEPTH));
        printf("4 %d %d %d %d\n", v-(HEIGHT*DEPTH+1), v-HEIGHT*DEPTH, v-(1+HEIGHT*DEPTH+DEPTH), v-1);
        printf("4 %d %d %d %d\n", v-(DEPTH+1), v-(1+HEIGHT*DEPTH+DEPTH), v-DEPTH, v-1); */
        printf("4 %d %d %d %d\n", v-(HEIGHT*DEPTH+DEPTH), v-HEIGHT*DEPTH, v-DEPTH, v-(1+HEIGHT*DEPTH+DEPTH));
        printf("4 %d %d %d %d\n", v-(HEIGHT*DEPTH+1), v-(1+HEIGHT*DEPTH+DEPTH), v-HEIGHT*DEPTH, v-1);
        printf("4 %d %d %d %d\n", v-(DEPTH+1), v-DEPTH, v-(1+HEIGHT*DEPTH+DEPTH), v-1);
        printf("4 %d %d %d %d\n", v-HEIGHT*DEPTH, v-DEPTH, v-1, v-(1+HEIGHT*DEPTH+DEPTH));
      }
  return 0;
}

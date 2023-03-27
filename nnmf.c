#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> // memset
#include <stdbool.h> // bool

// credits https://github.com/0xfffffff7/NMF
// compile with $ gcc nnmf.c -Wall -lm -o nnmf && ./nnmf

void matrix_add(double *c, const double *a, const double *b, const int n, const int m) {
  int i = 0;
  int j = 0;
  for (i = 0; i < n; i ++) {
    for (j = 0; j < m; j ++) {
      *(c + i * m + j) = *(a + i * m + j) + *(b + i * m + j);
    }
  }
}

int matrix_diff(const double *a, const double *b, const int n, const int m) {
  int i = 0;
  int j = 0;
  int diff = 0;
  for (i = 0; i < n; i ++) {
    for (j = 0; j < m; j ++) {
      diff = (int)pow(*(a + i * m + j) - *(b + i * m + j), 2);
    }
  }
  return diff;
}

void matrix_init(double *a, const int n, const int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      *(a + (i * m + j)) = (double)(rand() % 1000);
      srand((unsigned int)time(NULL) ^ (i + j));
    }
  }
}

void matrix_init_zero(double *a, const int n, const int m) {
  memset(a, 0x00, n * m * sizeof(double));
}

void matrix_copy(double *a, const double *b, const int n, const int m) {
  memcpy(a, b, n * m * sizeof(double));
}

void matrix_print(const double *a, const int n, const int m, int coefficient, bool bDisplayPoint) {
  if (bDisplayPoint) {
    for (int i = 0; i < n; i ++) {
      for (int j = 0; j < m; j ++) {
        printf("%12.1f ", *(a + (i * m + j)) * coefficient);
      }
      printf("\n");
    }
  }
  else{
    for (int i = 0; i < n; i ++) {
      for (int j = 0; j < m; j ++) {
        printf("%d ", (int)(*(a + (i * m + j)) * coefficient));
      }
      printf("\n");
    }
  }
  printf("\n");
}

void matrix_xl(double *c, const double *a, const double *b, const int l, const int m, const int n) {
  double t = 0;
  for (int i = 0; i < l; i ++) {
    for (int j = 0; j < n; j ++) {
      for (int k = 0; k < m; k ++) {
        t += *(a + (i * m + k)) * *(b + (k * n + j));
      }
      *(c + (i * n + j)) = t;
      t = 0;
    }
  }
}

void matrix_x(double *c, const double *a, const double *b, const int m, const int n) {
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      *(c + (i * n + j)) = *(a + (i * n + j)) * *(b + (i * n + j));
    }
  }
}

void matrix_divide(double *c, const double *a, const double *b, const int m, const int n) {
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      *(c + (i * n + j)) = *(a + (i * n + j)) / *(b + (i * n + j));
    }
  }
}

void matrix_transpose(const double *a, double *b, const int n, const int m) {
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < m; j ++) {
      *(b + j * n + i) = *(a + i * m + j);
    }
  }
}

void matrix_left_transpose_x(double *c, const double *a, const double *b, const int l, const int m, const int n) {
  double t[l * m];
  matrix_init_zero(t, l, m);
  matrix_transpose(a, t, l, m);
  matrix_xl(c, t, b, l, m, n);
}

void matrix_right_transpose_x(double *c, const double *a, const double *b, const int l, const int m, const int n) {
  double t[n * m];
  matrix_init_zero(t, n, m);
  matrix_transpose(b, t, n, m);
  matrix_xl(c, a, t, l, m, n);
}

void matrix_tx_left(double *c, const double *a, const double *b, const int l, const int m, const int n) {
  double t = 0;
  for (int i = 0; i < l; i ++) {
    for (int j = 0; j < n; j ++) {
      for (int k = 0; k < m; k ++) {
        t += *(a + (k * l + i)) * *(b + (k * n + j));
      }
      *(c + (i * n + j)) = t;
      t = 0;
    }
  }
}

void matrix_tx_right(double *c, const double *a, const double *b, const int l, const int m, const int n) {
  double t = 0;
  for (int i = 0; i < l; i ++) {
    for (int j = 0; j < n; j ++) {
      for (int k = 0; k < m; k ++) {
        t += *(a + (i * m + k)) * *(b + (j * m + k));
      }
      *(c + (i * n + j)) = t;
      t = 0;
    }
  }
}

void factorize(double *v, int row, int col, int features, unsigned int count) {
  double h[features * col];
  double hn[features * col];
  double hd[features * col];
  double hd_t[features * features];
  double w[row * features];
  double wn[row * features];
  double wd[row * features];
  double wd_t[row * col];
  double wh[row * col];

  matrix_init_zero((double*)wh, row, col);
  matrix_init(h, features, col);
  matrix_init_zero(hn, features, col);
  matrix_init_zero(hd, features, col);
  matrix_init_zero(hd_t, features, features);
  matrix_init(w, row, features);
  matrix_init_zero(wn, row, features);
  matrix_init_zero(wd, row, features);
  matrix_init_zero(wd_t, row, col);


  for (int i = 0; i < count; i ++) {
    matrix_xl(wh, w, h, row, features, col);
    int cost = matrix_diff(v, wh, row, col);

    if (cost == 0){
      break;
    }

    matrix_tx_left(hn, w, v, features, row, col);

    matrix_tx_left(hd_t, w, w, features, row, features);
    matrix_xl(hd, hd_t, h, features, features, col);

    matrix_x(h, h, hn, features, col);
    matrix_divide(h, h, hd, features, col);

    matrix_tx_right(wn, v, h, row, col, features);

    matrix_xl(wd_t, w, h, row, features, col);
    matrix_tx_right(wd, wd_t, h, row, col, features);

    matrix_x(w, w, wn, row, features);
    matrix_divide(w, w, wd, row, features);
  }


  printf("V =\n");
  matrix_print(v, row, col, 1, true);

  printf("W x H =\n");
  matrix_print(wh, row, col, 1, true);

  printf("W =\n");
  matrix_print(w, row, features, 1, true);

  printf("H * 1000 =\n");
  matrix_print(h, features, col, 1000, true);
}


int main(int argc, char* argv[]) {
  int FEATURES = 2;
  int ROW = 5;
  int COL = 4;
  int COUNT = 50000;

  const double nARRAY[5][4] = {
    {  50.0,  50.0,  50.0,  50.0 },
    { 100.0,  20.0,  30.0,  90.0 },
    { 100.0,  50.0,  40.0, 100.0 },
    {  90.0,  10.0,  10.0,  80.0 },
    {  10.0, 100.0,  90.0,  20.0 }
  };

  factorize((double*)nARRAY, ROW, COL, FEATURES, COUNT);

  return 0;
}


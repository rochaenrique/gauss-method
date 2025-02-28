#include "matrix.h"

int main()
{
  Matrix *mat = matrix_alloc(5, 5);
  matrix_randomize(mat);
  printf("det() = %g\n", matrix_gauss_det(mat));
  matrix_free(mat);
  
  return 0;
}

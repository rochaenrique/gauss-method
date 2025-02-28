#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define ROUND_MAX 0.0001f

typedef struct {
  float **buf;
  size_t m, n;
} Matrix;

#define matrix_foreach(mat, func)				\
  do {											\
	for (size_t i = 0; i < mat->m; i++)			\
	  for (size_t j = 0; j < mat->n; j++)		\
		(func)(&mat->buf[i][j]);				\
  } while (0);

#define matrix_foreach_set(mat, func)			\
  do {											\
	for (size_t i = 0; i < mat->m; i++)			\
	  for (size_t j = 0; j < mat->n; j++)		\
		mat->buf[i][j] = (func)();				\
  } while (0);									

#define matrix_set(mat, value)				    \
  do {											\
	for (size_t i = 0; i < (mat)->m; i++)		\
	  for (size_t j = 0; j < (mat)->n; j++)		\
		(mat)->buf[i][j] = (value);				\
  } while (0);									

#define matrix_zero(mat) matrix_set(mat, 0)
#define rand_normalized(max, min) (float)(min + (max - min) * ((float)rand() / (float)RAND_MAX))
#define _is_round_zero(num) ((num) >= -ROUND_MAX && (num) <= ROUND_MAX)

// declarations
Matrix *matrix_alloc(size_t, size_t);
void matrix_free(Matrix *);
void matrix_randomize(Matrix *);
void matrix_gauss(Matrix *);
void submatrix_gauss(Matrix *, size_t, size_t);
void matrix_gauss_jordan(Matrix *);
float matrix_gauss_det(Matrix *);
size_t matrix_find_not_zero(Matrix *mat, size_t col, size_t from);
void matrix_swap_lines(Matrix *, size_t, size_t);
void matrix_print(Matrix *);
void matrix_round_zero(Matrix *);

// implementations
Matrix *matrix_alloc(size_t m, size_t n)
{
  Matrix *matrix = malloc(sizeof(Matrix));
  matrix->buf = calloc(m, sizeof(int *));
  matrix->m = m;
  matrix->n = n;
  for (; m > 0; --m) 
	matrix->buf[m-1] = calloc(n, sizeof(int));
  
  return matrix;
}

void matrix_free(Matrix *mat)
{
  for (size_t i = 0; i < mat->m; i++) 
	free(mat->buf[i]);
  
  free(mat->buf);
  mat->m = 0;
  mat->n = 0;
  free(mat);
}

int _rand() { return rand_normalized(-20, 20); }
void matrix_randomize(Matrix *mat)
{
  matrix_foreach_set(mat, _rand);
}

void matrix_gauss(Matrix *mat)
{
  submatrix_gauss(mat, 0, 0);
}

void submatrix_gauss(Matrix *mat, size_t row, size_t col)
{
  if (!(row < mat->m && col < mat->n))
	return;
  
  // si hay un 0, intercambiamos
  if (mat->buf[row][col] == 0.0f) { 
	size_t swap = matrix_find_not_zero(mat, row, col);
	if (swap != 0)
	  matrix_swap_lines(mat, row, swap);
  }
  
  // hacer ceros
  float alpha; 
  for (size_t i = row+1; i < mat->m; i++) {
	alpha = -mat->buf[i][col] / mat->buf[row][col];
	
	for (size_t j = col; j < mat->n; j++) 
	  mat->buf[i][j] += alpha * mat->buf[row][j];
  }

  submatrix_gauss(mat, row+1, col+1);
}

void matrix_gauss_jordan(Matrix *mat)
{
  matrix_gauss(mat);
  
  // hacer unos
  int row;
  size_t col;
  float alpha;
  for (row = mat->m - 1; row > -1; row--) {
	for (col = 0; col < mat->n && _is_round_zero(mat->buf[row][col]); col++);
	alpha = 1.0f / mat->buf[row][col];

	for (size_t i = col; i < mat->n; i++) 
	  mat->buf[row][i] *= alpha;

	// hacer ceros arriba del uno
	for (int i = row - 1; i > -1; i--) {
	  alpha = -mat->buf[i][col] / mat->buf[row][col];
	  for (size_t j = col; j < mat->n; j++) 
		mat->buf[i][j] += alpha * mat->buf[row][j];
	}
  }
}

float matrix_gauss_det(Matrix *mat)
{
  assert(mat->m == mat->n && "Only square matrices have determinants");
  matrix_gauss(mat);

  float resul = 1.0f;
  for (size_t i = 0; i < mat->m; i++) 
	resul *= mat->buf[i][i];

  return resul;
}

size_t matrix_find_not_zero(Matrix *mat, size_t col, size_t from)
{
  size_t toswap = from;
  while (toswap < mat->m && mat->buf[toswap][col] == 0)
	toswap++;
  
  return (mat->buf[toswap][col] != 0) ? toswap : 0;
}

void matrix_swap_lines(Matrix *mat, size_t a, size_t b)
{
  assert(a < mat->m && b < mat->m && "Lines to swap are not in bounds");
  float *tmp = mat->buf[a];
  mat->buf[a] = mat->buf[b];
  mat->buf[b] = tmp;
}

void matrix_print(Matrix *mat)
{
  printf("Matrix %zux%zu\n", mat->m, mat->n);
  for (size_t i = 0; i < mat->m; i++) {
	printf("| ");
	for (size_t j = 0; j < mat->n; j++)
	  printf("%8g ", mat->buf[i][j]);
	printf(" | \n");
  }
}

void _round_if(float *a)
{
  if (_is_round_zero(*a))
	*a = 0.0f;
}

void matrix_round_zero(Matrix *mat)
{
  matrix_foreach(mat, _round_if);
}


#endif // MATRIX_H_

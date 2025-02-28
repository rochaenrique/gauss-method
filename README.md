# Gauss Method
A simple implementation of the gaussian-method to solve systems of linear equations and calculate determinants.
The file `matrix.h` consists of the library in, both portable to nativ and wasm executables.

## Usage
Check `/examples`

`/examples/basic.c`
```c
  Matrix *mat = matrix_alloc(5, 5);
  matrix_randomize(mat);
  printf("det() = %g\n", matrix_gauss_det(mat));
  matrix_free(mat);
```

## Compiling to Wasm
The build.sh script compiles the basic.c example nativly and to wasm (https://surma.dev/things/c-to-webassembly/);

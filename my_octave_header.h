///Copyright 2021 311CA Vlad Negoita <vlad1negoita@gmail.com>

#ifndef MY_OCTAVE_HEADER_H
#define MY_OCTAVE_HEADER_H

///**a = the address to the elements of the matrix
///lines = the number of lines
///columns = the number of columns
///sum = the sum of the elements of a matrix
typedef struct {
	int **a;
	int lines, columns;
	int sum;
} matrix;

int modulo_add(int x, int y);
matrix alloc_matrix(int n, int m);
void read_matrix(matrix *v, int index);
void merge_sort(matrix *v, int left, int right);
void free_matrix(matrix mat);
void transpose(matrix *v, int index);
void resize(matrix *v, int matrix_cnt);
matrix multiply(matrix *v, int index1, int index2);
matrix matrix_add(matrix A, matrix B, int substract);
void matrix_to_blocks(matrix A, matrix *A1);
matrix strassen(matrix A, matrix B);
void print_matrix(matrix *v, int index);
void increase_size(matrix **v, int *v_size, int *matrix_cnt, int dim);
void decrease_size(matrix **v, int *v_size, int *matrix_cnt, int dim);

#endif

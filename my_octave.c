///Copyright 2021 311CA Vlad Negoita <vlad1negoita@gmail.com>
#include <stdio.h>
#include <stdlib.h>
#include "my_octave_header.h"

#define MOD 10007
#define INDEX_ERROR "No matrix with the given index"
#define COMMAND_ERROR "Unrecognized command"
#define MULTIPLICATION_ERROR "Cannot perform matrix multiplication"

int main(void)
{
	int matrix_cnt = 0, v_size = 0, finished = 0;
	int dim = sizeof(matrix);
	matrix *v = NULL;
	char operation;
	while (!finished) {
		scanf("%c", &operation);///next operation
		if (operation == 'L') {///new matrix
			increase_size(&v, &v_size, &matrix_cnt, dim);
			int n, m;
			scanf("%d %d", &n, &m);
			v[matrix_cnt - 1] = alloc_matrix(n, m);
			read_matrix(v, matrix_cnt - 1);
		} else if (operation == 'D') {///dimensions
			int index;
			scanf("%d", &index);
			if (index >= matrix_cnt || index < 0)
				printf("%s\n", INDEX_ERROR);
			else
				printf("%d %d\n", v[index].lines, v[index].columns);
		} else if (operation == 'P') {///print
			int index;
			scanf("%d", &index);
			if (index >= matrix_cnt || index < 0)
				printf("%s\n", INDEX_ERROR);
			else
				print_matrix(v, index);
		} else if (operation == 'C') {///cut / resize
			resize(v, matrix_cnt);
		} else if (operation == 'M' || operation == 'S') {
			int index1, index2;
			scanf("%d %d", &index1, &index2);
			if (index1 >= matrix_cnt || index1 < 0) {
				printf("%s\n", INDEX_ERROR);
			} else if (index2 >= matrix_cnt || index2 < 0) {
				printf("%s\n", INDEX_ERROR);
			} else if (v[index1].columns != v[index2].lines) {
				printf("%s\n", MULTIPLICATION_ERROR);
			} else if (operation == 'M') {
				increase_size(&v, &v_size, &matrix_cnt, dim);
				v[matrix_cnt - 1] = multiply(v, index1, index2);
			} else {
				increase_size(&v, &v_size, &matrix_cnt, dim);
				v[matrix_cnt - 1] = strassen(v[index1], v[index2]);
			}
		} else if (operation == 'O') {///sort
			merge_sort(v, 0, matrix_cnt - 1);
		} else if (operation == 'T') {///transpose
			int index;
			scanf("%d", &index);
			if (index >= matrix_cnt)
				printf("%s\n", INDEX_ERROR);
			else
				transpose(v, index);
		} else if (operation == 'F') {///delete a matrix
			int index;
			scanf("%d", &index);
			if (index >= matrix_cnt || index < 0) {
				printf("%s\n", INDEX_ERROR);
			} else {
				free_matrix(v[index]);
				for (int i = index; i < matrix_cnt - 1; ++i)
					v[i] = v[i + 1];
				decrease_size(&v, &v_size, &matrix_cnt, dim);
			}
		} else if (operation == 'Q') {///quit
			for (int index = 0; index < matrix_cnt; ++index)
				free_matrix(v[index]);
			matrix_cnt = 0;
			v_size = 0;
			free(v);
			finished = 1;///exit
		} else {
			printf("%s\n", COMMAND_ERROR);
		}
		scanf("%c", &operation);///residual '\n'
	}
	return 0;
}

///the correct reminder of (x + y) % MOD
int modulo_add(int x, int y)
{
	return (((x + y) % MOD + MOD) % MOD);
}

///alocates memory for a matrix
matrix alloc_matrix(int n, int m)
{
	matrix new_matrix;
	new_matrix.lines = n;
	new_matrix.columns = m;
	new_matrix.a = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; ++i)
		new_matrix.a[i] = (int *)malloc(m * sizeof(int));
	new_matrix.sum = 0;
	return new_matrix;
}

///read values in a matrix
void read_matrix(matrix *v, int index)
{
	for (int i = 0; i < v[index].lines; ++i) {
		for (int j = 0; j < v[index].columns; ++j) {
			scanf("%d", &v[index].a[i][j]);
			v[index].sum = modulo_add(v[index].sum, v[index].a[i][j]);
		}
	}
}

///sort the matrix list using merge_sort algorithm
void merge_sort(matrix *v, int left, int right)
{
	if (left >= right)
		return;
	int middle = (left + right) / 2;
	merge_sort(v, left, middle);
	merge_sort(v, middle + 1, right);
	matrix *aux;
	aux = (matrix *)malloc((right - left + 1) * sizeof(matrix));
	int i = left, j = middle + 1, k = 0;
	while (i <= middle && j <= right) {
		if (v[i].sum < v[j].sum)
			aux[k++] = v[i++];
		else
			aux[k++] = v[j++];
	}
	while (i <= middle)
		aux[k++] = v[i++];
	while (j <= right)
		aux[k++] = v[j++];
	for (int ind = left; ind <= right; ++ind)
		v[ind] = aux[ind - left];
	free(aux);
}

///free all the dynamic allocated memory from a matrix
void free_matrix(matrix mat)
{
	for (int i = 0; i < mat.lines; ++i)
		free(mat.a[i]);
	free(mat.a);
}

///transpose a matrix
void transpose(matrix *v, int index)
{
	matrix aux = alloc_matrix(v[index].columns, v[index].lines);
	aux.sum = v[index].sum;
	for (int i = 0; i < v[index].lines; ++i)
		for (int j = 0; j < v[index].columns; ++j)
			aux.a[j][i] = v[index].a[i][j];
	free_matrix(v[index]);
	v[index] = aux;
}

/// cut the matrix so that the elements are
///at the intersection of lines[i] and columns[j]
void resize(matrix *v, int matrix_cnt)
{
	int index;
	scanf("%d", &index);
	int n, m;
	int *lines, *columns;
	scanf("%d", &n);
	lines = (int *)malloc(n * sizeof(int));
	for (int i = 0; i < n; ++i)
		scanf("%d", &lines[i]);
	scanf("%d", &m);
	columns = (int *)malloc(m * sizeof(int));
	for (int j = 0; j < m; ++j)
		scanf("%d", &columns[j]);
	if (index < 0 || index >= matrix_cnt) {
		printf("%s\n", INDEX_ERROR);
		free(lines);
		free(columns);
		return;
	}
	matrix resized = alloc_matrix(n, m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			resized.a[i][j] = v[index].a[lines[i]][columns[j]];
			resized.sum = modulo_add(resized.sum, resized.a[i][j]);
		}
	}
	free_matrix(v[index]);
	free(lines);
	free(columns);
	v[index] = resized;
}

///returns v[index1] * v[index2] (matrix product)
matrix multiply(matrix *v, int index1, int index2)
{
	matrix aux;
	aux = alloc_matrix(v[index1].lines, v[index2].columns);
	for (int i = 0; i < v[index1].lines; ++i)
		for (int j = 0; j < v[index2].columns; ++j) {
			int s = 0;
			for (int k = 0; k < v[index1].columns; ++k)
				s = modulo_add(s, v[index1].a[i][k] * v[index2].a[k][j] % MOD);
			aux.a[i][j] = s;
			aux.sum = modulo_add(aux.sum, s);
		}
	return aux;
}

///substract = 0 => returns A + B
///substract = 1 => returns A - B
matrix matrix_add(matrix A, matrix B, int substract)
{
	matrix aux;
	aux = alloc_matrix(A.lines, B.columns);
	for (int i = 0; i < aux.lines; ++i) {
		for (int j = 0; j < aux.columns; ++j) {
			if (!substract)
				aux.a[i][j] = modulo_add(A.a[i][j], B.a[i][j]);
			else
				aux.a[i][j] = modulo_add(A.a[i][j], -B.a[i][j]);
		}
	}
	aux.sum = modulo_add(A.sum, B.sum);
	return aux;
}

///useful for strassen
///copy the values from A to A1 ... A4
void matrix_to_blocks(matrix A, matrix *A1)
{
	int dim = A.lines / 2;
	for (int i = 0; i < A.lines; ++i)
		for (int j = 0; j < A.lines; ++j)
			if (i < dim && j < dim)
				A1[0].a[i][j] = A.a[i][j];
			else if (i < dim)
				A1[1].a[i][j - dim] = A.a[i][j];
			else if (j < dim)
				A1[2].a[i - dim][j] = A.a[i][j];
			else
				A1[3].a[i - dim][j - dim] = A.a[i][j];
}

///returns A * B
///works only if A and B have the dimensions a power of 2
matrix strassen(matrix A, matrix B)
{
	int dim = A.lines / 2;
	matrix aux, *A1, *B1, *M1;
	if (A.lines == 1) {///exit recursion condition
		aux = alloc_matrix(1, 1);
		aux.a[0][0] =  A.a[0][0] * B.a[0][0];
		aux.sum = A.a[0][0] * B.a[0][0];
		return aux;
	}
	A1 = (matrix *)malloc(4 * sizeof(matrix));
	B1 = (matrix *)malloc(4 * sizeof(matrix));
	for (int i = 0; i < 4; ++i) {
		A1[i] = alloc_matrix(dim, dim);
		B1[i] = alloc_matrix(dim, dim);
	}
	matrix_to_blocks(A, A1);///the blocks A1, A2, A3, A4
	matrix_to_blocks(B, B1);///the blocks B1, B2, B3, B4
	matrix tmp1 = matrix_add(A1[0], A1[3], 0);
	matrix tmp2 = matrix_add(B1[0], B1[3], 0);
	M1 = (matrix *)malloc(7 * sizeof(matrix));
	M1[0] = strassen(tmp1, tmp2);
	free_matrix(tmp1);
	free_matrix(tmp2);
	M1[1] = strassen(tmp1 = matrix_add(A1[2], A1[3], 0), B1[0]);
	free_matrix(tmp1);
	M1[2] = strassen(A1[0], tmp1 = matrix_add(B1[1], B1[3], 1));
	free_matrix(tmp1);
	M1[3] = strassen(A1[3], tmp1 = matrix_add(B1[2], B1[0], 1));
	free_matrix(tmp1);
	M1[4] = strassen(tmp1 = matrix_add(A1[0], A1[1], 0), B1[3]);
	free_matrix(tmp1);
	tmp1 = matrix_add(A1[2], A1[0], 1);
	tmp2 = matrix_add(B1[0], B1[1], 0);
	M1[5] = strassen(tmp1, tmp2);
	free_matrix(tmp1);
	free_matrix(tmp2);
	tmp1 = matrix_add(A1[1], A1[3], 1);
	tmp2 = matrix_add(B1[2], B1[3], 0);
	M1[6] = strassen(tmp1, tmp2);
	free_matrix(tmp1);
	free_matrix(tmp2);
	for (int i = 0; i < 4; ++i) {
		free_matrix(A1[i]);
		free_matrix(B1[i]);
	}
	free(A1);
	free(B1);
	aux = alloc_matrix(A.lines, A.columns);
	for (int i = 0; i < aux.lines; ++i)
		for (int j = 0; j < aux.columns; ++j) {
			int k = i - dim, l = j - dim;
			if (i < dim && j < dim) {
				aux.a[i][j] = modulo_add(M1[0].a[i][j], M1[3].a[i][j]);
				aux.a[i][j] = modulo_add(aux.a[i][j], -M1[4].a[i][j]);
				aux.a[i][j] = modulo_add(aux.a[i][j], M1[6].a[i][j]);
			} else if (i < dim) {
				aux.a[i][j] = modulo_add(M1[2].a[i][l], M1[4].a[i][l]);
			} else if (j < dim) {
				aux.a[i][j] = modulo_add(M1[1].a[k][j], M1[3].a[k][j]);
			} else {
				aux.a[i][j] = modulo_add(M1[0].a[k][l], -M1[1].a[k][l]);
				aux.a[i][j] = modulo_add(aux.a[i][j], M1[2].a[k][l]);
				aux.a[i][j] = modulo_add(aux.a[i][j], M1[5].a[k][l]);
			}
			aux.sum = modulo_add(aux.sum, aux.a[i][j]);
		}
	for (int i = 0; i < 7; ++i)
		free_matrix(M1[i]);
	free(M1);
	return aux;
}

void print_matrix(matrix *v, int index)
{
	for (int i = 0; i < v[index].lines; ++i) {
		for (int j = 0; j < v[index].columns; ++j)
			printf("%d ", v[index].a[i][j]);
		printf("\n");
	}
}

///resizable array -> the length is a power of 2 (2 ^ k)
///this function increases the size (doubles it) when needed
///also incrementing the matrix_cnt
void increase_size(matrix **v, int *v_size, int *matrix_cnt, int dim)
{
	*matrix_cnt = *matrix_cnt + 1;
	if (!(*v)) {
		*v = (matrix *)malloc(*matrix_cnt * dim);
		*v_size = 1;
	} else if (*v_size < *matrix_cnt) {
		*v_size = *v_size * 2;
		*v = (matrix *)realloc(*v, *v_size * dim);
	}
}

///resizable array -> the length is a power of 2 (2 ^ k)
///this function decreases the size (halves it) when possible
///also decrementing the matrix_cnt
void decrease_size(matrix **v, int *v_size, int *matrix_cnt, int dim)
{
	*matrix_cnt = *matrix_cnt - 1;
	if (matrix_cnt == 0) {
		free(*v);
	} else if (*v_size > *matrix_cnt * 2) {
		*v_size = *v_size / 2;
		*v = (matrix *)realloc(*v, *v_size * dim);
	}
}

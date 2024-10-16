#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CORRECT 0
/* Ошибка, некорректная матрица */
#define INCOR_MATR 1
/* Ошибка вычисления (несовпадающие размеры матриц; матрица, для которой нельзя
   провести вычисления и т.д.) */
#define CALC_ER 2

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

void Print_Matrix(matrix_t A);
int is_True(matrix_t *matrix);

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_eq_matrix(matrix_t *A, matrix_t *B);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
double s21_determinant_recursive(matrix_t *A);
double s21_get_determinant(matrix_t *A);

int Minor(int row, int column, matrix_t *A, matrix_t *result);
int calc_helper(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

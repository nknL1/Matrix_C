#include "s21_matrix.h"

// 0 - OK
// 1 - Ошибка, некорректная матрица
// 2 - Ошибка вычисления (несовпадающие размеры матриц;
//  матрица, для которой нельзя провести вычисления и т.д.)

// void Print_Matrix(matrix_t A) {
//   for (int i = 0; i < A.rows; i++) {
//     for (int j = 0; j < A.columns; j++) {
//       printf("%f ", A.matrix[i][j]);
//     }
//     printf("\n");
//   }
// }
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = CORRECT;
  if (rows < 1 || columns < 1) {
    error = INCOR_MATR;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    error = INCOR_MATR;
    if (result->matrix != NULL) {
      error = CORRECT;
      for (int i = 0; i < result->rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
        if (result->matrix[i] == NULL) {
          error = INCOR_MATR;  // Ошибка, некорректная матрица
        }
      }
    }
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i]) {
        free(A->matrix[i]);
      }
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A->rows) {
    A->rows = 0;
  }
  if (A->columns) {
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error = FAILURE;
  if (is_True(A) && is_True(B)) {
    error = FAILURE;
    if (A->rows == B->rows && A->columns == B->columns) {
      error = SUCCESS;
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-6) {
            error = FAILURE;
          }
        }
      }
    }
  }
  return error;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A) && is_True(B)) {
    error = CALC_ER;
    if (A->rows == B->rows && A->columns == B->columns) {
      error = s21_create_matrix(A->rows, A->columns, result);
      if (error == CORRECT) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
          }
        }
      }
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A) && is_True(B)) {
    error = CALC_ER;
    if (A->rows == B->rows && A->columns == B->columns) {
      error = s21_create_matrix(A->rows, A->columns, result);
      if (error == CORRECT) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
          }
        }
      }
    }
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A)) {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == CORRECT) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A) && is_True(B)) {
    if (A->columns == B->rows) {
      error = s21_create_matrix(A->rows, B->columns, result);
      if (error == CORRECT) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < B->columns; j++) {
            for (int k = 0; k < B->rows; k++) {
              result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
            }
          }
        }
      }
    } else {
      error = CALC_ER;
    }
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = INCOR_MATR;
  if (A->columns && A->rows) {
    error = s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return error;
}

int is_True(matrix_t *matrix) {
  int res;
  if (matrix != NULL && matrix->matrix != NULL && matrix->rows >= 1 &&
      matrix->columns >= 1) {
    res = 1;
  } else {
    res = 0;
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = INCOR_MATR;
  if (is_True(A)) {
    error = CALC_ER;
    if (A->rows == A->columns) {
      error = CORRECT;
      *result = A->matrix[0][0];
      if (A->rows != 1) {
        *result = s21_determinant_recursive(A);
      }
    }
  }
  return error;
}

double s21_determinant_recursive(matrix_t *A) {
  double result = 0;
  if (A->rows == 2) {
    result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    for (int i = 0; i < A->rows; i++) {
      matrix_t minor;
      Minor(1, i + 1, A, &minor);
      result +=
          pow((-1), i) * A->matrix[0][i] * s21_determinant_recursive(&minor);
      s21_remove_matrix(&minor);
    }
  }
  return result;
}

int Minor(int row, int column, matrix_t *A, matrix_t *result) {
  int error = 1;
  if (A->matrix != NULL) {
    error = s21_create_matrix(A->rows - 1, A->columns - 1, result);
    if (error == 0) {
      int m, n;
      for (int i = 0; i < A->rows; i++) {
        m = i;
        // Смещение индекса строки m, если текущая строка совпадает с заданной
        // row
        if (i > row - 1) {
          m--;
        }
        // Смещение индекса столбца n, если текущий столбец совпадает с заданным
        // column
        for (int j = 0; j < A->columns; j++) {
          n = j;
          if (j > column - 1) {
            n--;
          }
          // Копирую элемента из A в result, исключая строку row и столбец
          // column
          if (i != row - 1 && j != column - 1) {
            result->matrix[m][n] = A->matrix[i][j];
          }
        }
      }
    }
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A)) {
    error = CALC_ER;
    if (A->rows == A->columns) {
      error = s21_create_matrix(A->rows, A->columns, result);
      if (error == CORRECT) {
        error = calc_helper(A, result);
      }
    }
  }
  return error;
}

int calc_helper(matrix_t *A, matrix_t *result) {
  int error = 0;
  result->matrix[0][0] = 1;
  if (A->rows != 1) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        double deted;
        matrix_t minored;
        error = Minor(i + 1, j + 1, A, &minored);
        if (error == 0) {
          error = s21_determinant(&minored, &deted);
          if (error == 0) {
            result->matrix[i][j] = pow((-1), i + j) * deted;
          }
        }
        s21_remove_matrix(&minored);
      }
    }
  } else {
    error = 2;
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = INCOR_MATR;
  if (is_True(A)) {
    error = CALC_ER;
    double det;
    s21_determinant(A, &det);
    if (fabs(det - 0) > 1e-6) {
      matrix_t tmp_calc;
      error = s21_calc_complements(A, &tmp_calc);
      if (error == CORRECT) {
        matrix_t tmp_trans;
        error = s21_transpose(&tmp_calc, &tmp_trans);
        if (error == CORRECT) {
          s21_mult_number(&tmp_trans, 1 / det, result);
        }
        s21_remove_matrix(&tmp_trans);
      }
      s21_remove_matrix(&tmp_calc);
    }
  }
  return error;
}
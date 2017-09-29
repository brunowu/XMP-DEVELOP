#include "vector.h"
#include "matrix.h"

#ifndef MVECMATRIX_H
#define MVECMATRIX_H

#define VECTOR_DIVID(vec, x) vector_divid(&vec, x)
#define VECTOR_MULTIPLE(vec, x) vector_multiple(&vec, x)
#define VECTOR_ADD_VECTOR(vec1, vec2, newV) vector_add_vector(&vec1, &vec2, &newV);
#define VECTOR_REDUCE_VECTOR(vec1, vec2, newV) vector_reduce_vector(&vec1, &vec2, &newV)
#define VECTOR_INNER_PRODUIT(vec1, vec2, h) vector_inner_produit(&vec1, &vec2, &h)
#define VECTOR_ABS(vec, h) vector_abs(&vec, &h)
#define STOCK_VECTOR_IN_MATRIX(mat, vec, ligne, colone) stock_vector_in_matrix(&mat, &vec, ligne, colone)
#define STOCK_MVECTOR_IN_MATRIX(mat, v, ligne) stock_mvector_in_matrix(&mat, &&v, ligne)
#define MATRIX_GET_VECTOR(mat, vec, ligne) matrix_get(&mat, &vec, ligne)
#define MATRIX_MULTIPLE_VECTOR(mat, vec, newV) matrix_multiple_vector(&mat, &vec, &newV)
#define MATRIX_MULTIPLE_MATRIX(mat1, mat2, newMat) matrix_multiple_matrix(&mat1, &mat2, &newMat)
#define UPPER_TRIANGLE_MATRIX_INVERSE(mat, mati) upper_triangle_matrix_inverse(&mat, &mati)

void vector_divid(vector * v, double x);
void vector_multiple(vector * v, double x);
void vector_add_vector(vector * v1, vector * v2, vector * newV);
void vector_reduce_vector(vector * v1, vector * v2, vector * newV);
void vector_inner_produit(vector * v1, vector * v2, double * h);
void vector_abs(vector * v, double * h);
void stock_vector_in_matrix(matrix * mat, vector * v, int ligne);
void stock_mvector_in_matrix(matrix * mat, void ** v, int ligne, int colone);
void matrix_get_vector(matrix * mat, vector * v, int ligne);
void matrix_multiple_vector(matrix * mat, vector * v, vector * newV);
void matrix_multiple_matrix(matrix * mat1, matrix * mat2, matrix * newMat);
void upper_triangle_matrix_inverse(matrix * mat, matrix * mati);
 
/********************
Matrix Ellpack
*********************/
#define MATRIX_ELL_MULTIPLE_VECTOR(mat, vec, newV) matrix_multiple_vector(&mat, &vec, &newV)
#define MATRIX_ELL_MULTIPLE_MATRIX(mat1, mat2, newMat) matrix_multiple_matrix(&mat1, &mat2, &newMat)
#define MATRIX_CONVERT_MATRIX_ELL(mat, matell) matrix_convert_matrix_ell(&mat, &matell)
void matrix_ell_multiple_vector(matrix_ELL * mat, vector * v, vector * newV);
void matrix_ell_multiple_matrix(matrix_ELL * mat1, matrix * mat2, matrix * newMat);
void matrix_convert_matrix_ell(matrix * mat, matrix_ELL * matell);
#endif
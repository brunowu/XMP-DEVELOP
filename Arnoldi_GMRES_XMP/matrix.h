#include "vector.h"

#ifndef MATRIX_H
#define MATRIX_H

#define MATRIX_INIT(mat, ligne, colone) matrix_init(&mat, ligne, colone)
#define MATRIX_INIT_TOTAL(mat, start, end) matrix_init_total(&mat, start, end)
#define MATRIX_COMPLETE_LIGNE(mat) matrix_complete_ligne(&mat)
#define MATRIX_RESIZE(mat, ligne, colone) matrix_resize(&mat, ligne, colone)
#define MATRIX_RESIZE_LIGNE(mat, ligne) matrix_resize_ligne(&mat, ligne)
#define MATRIX_RESIZE_COLONE(mat, colone) matrix_resize_colone(&mat, colone)
#define MATRIX_ADD(mat, ligne, item) matrix_add(&mat, ligne, (void *) item)
#define MATRIX_ADD_DUPLICATE(mat, ligne, item) matrix_add_duplicate(&mat, ligne, (void *) item)
#define MATRIX_COPY(mat1, mat2) matrix_copy(&mat1, &mat2)
#define MATRIX_DUPLICATE(mat1, mat2) matrix_duplicate(&mat1, &mat2)
#define MATRIX_SET(mat, ligne, colone, item) matrix_set(&mat, ligne, colone, (void *) item)
#define MATRIX_TRANSPOSE(mat1, mat2) matrix_transpose(&mat1, &mat2);
#define MATRIX_GET(mat, type, ligne, colone) matrix_get(&mat, ligne, colone)
#define MATRI_GET_LIGNE(mat) matrix_get_ligne(&mat)
#define MATRI_GET_CONLONE(mat) matrix_get_colone(&mat)
#define MATRIX_TOTAL_LIGNE(mat, ligne) matrix_total_ligne(&mat, ligne)
#define MATRIX_TOTAL(mat) matrix_total(&mat);
#define MATRIX_ITEMS(mat) matrix_items(&mat);
#define MATRIX_DELETE(mat, ligne, colone) matrix_delete(&mat, ligne, colone)
#define MATRIX_DELETE_LIGNE(mat, newm, ligne) matrix_delete_ligne(&mat, &newm, ligne)
#define MATRIX_DELETE_COLONE(mat, newm, colone) matrix_delete_colone(&mat, &newm, colone)
#define MATRIX_FREE(mat) matrix_free(&mat)
#define MATRIX_SHOW(mat) matrix_show(&mat)

typedef struct matrix {
	void *** items;
	int ligne;
	int colone;
	int * total;
} matrix;

void matrix_init(matrix * mat, int ligne, int colone);
static void matrix_init_total(matrix * mat, int start, int end);
void matrix_complete_ligne(matrix * mat);
int matrix_total_ligne(matrix * mat, int ligne);
void matrix_resize(matrix * mat, int ligne, int colone);
void matrix_resize_ligne(matrix * mat, int ligne);
void matrix_resize_colone(matrix * mat, int colone);
int matrix_get_ligne(matrix * mat);
int matrix_get_colone(matrix * mat);
int * matrix_total(matrix * mat);
void *** matrix_items(matrix * mat);
void matrix_add(matrix * mat, int ligne, void * item);
void matrix_add_duplicate(matrix * mat, int ligne, void * item);
void matrix_copy(matrix * mat1, matrix * mat2);
void matrix_duplicate(matrix * mat1, matrix * mat2);
void matrix_set(matrix * mat, int ligne, int colone, void * item);
void matrix_transpose(matrix * mat1, matrix * mat2);
void * matrix_get(matrix * mat, int ligne, int colone);
void matrix_delete(matrix * mat, int ligne, int colone);
void matrix_delete_ligne(matrix * mat, matrix * newm, int ligne);
void matrix_delete_colone(matrix * mat, matrix * newm, int colone);
void matrix_free(matrix * mat);
void matrix_show(matrix * mat);

/*******************************
Matrix Ellpack
********************************/
#define MATRIX_ELL_INIT(mat, ligne, colone, n_colone) matrix_init(&mat, ligne, colone, n_colone)
#define MATRIX_ELL_INIT_TOTAL(mat, start, end) matrix_init_total(&mat, start, end)
#define MATRIX_ELLCOMPLETE_LIGNE(mat) matrix_complete_ligne(&mat)
#define MATRIX_ELL_TOTAL_LIGNE(mat, ligne) matrix_total_ligne(&mat, ligne)
#define MATRIX_ELL_GET_LIGNE(mat) matrix_get_ligne(&mat)
#define MATRIX_ELL_GET_CONLONE(mat) matrix_get_colone(&mat)
#define MATRIX_ELL_N_COLONE(mat) matrix_ell_n_colone(&mat)
#define MATRIX_ELL_TOTAL(mat) matrix_ell_total(&mat)
#define MATRIX_ELL_ITEMS(mat) matrix_ell_items(&mat)
#define MATRIX_ELL_COPY(mat1, mat2) matrix_copy(&mat1, &mat2)
#define MATRIX_ELL_DUPLICATE(mat1, mat2) matrix_duplicate(&mat1, &mat2)
#define MATRIX_ELL_SET(mat, ligne, colone, item) matrix_set(&mat, ligne, colone, (void *) item)
#define MATRIX_ELL_GET_ITEMS(mat) matrix_ell_get_items(&mat)
#define MATRIX_ELL_GET(mat, ligne, colone) (type) matrix_get(&mat, ligne, colone)
#define MATRIX_ELL_FREE(mat) matrix_free(&mat)
#define MATRIX_ELL_SHOW(mat) matrix_show(&mat)
#define MATRIX_ELL_ADD(mat, ligne, colone, item) matrix_ell_add(&mat, ligne, colone, (void *) item)
#define MATRIX_ELL_ADD_DUPLICATE(mat, ligne, colone, item) matrix_ell_add_duplicate(&mat, ligne, colone, (void *) item)

typedef struct matrix_ELL {
	void *** items;
	int n_colone;
	int ligne;
	int colone;
	int * total;
} matrix_ELL;

void matrix_ell_init(matrix_ELL * mat, int ligne, int colone, int n_colone);
static void matrix_ell_init_total(matrix_ELL * mat, int start, int end);
void matrix_ell_complete_ligne(matrix_ELL * mat);
int matrix_ell_total_ligne(matrix_ELL * mat, int ligne);
int matrix_ell_get_ligne(matrix_ELL * mat);
int matrix_ell_get_colone(matrix_ELL * mat);
int matrix_ell_n_colone(matrix_ELL * mat);
int * matrix_ell_total(matrix_ELL * mat);
void *** matrix_ell_items(matrix_ELL * mat);
void matrix_ell_add(matrix_ELL * mat, int ligne, int colone, void * item);
void matrix_ell_add_duplicate(matrix_ELL * mat, int ligne, int colone, void * item);
void matrix_ell_copy(matrix_ELL * mat1, matrix_ELL * mat2);
void matrix_ell_duplicate(matrix_ELL * mat1, matrix_ELL * mat2);
void matrix_ell_set(matrix_ELL * mat, int ligne, int colone, void * item);
void * matrix_ell_get(matrix_ELL * mat, int ligne, int colone);

void matrix_ell_free(matrix_ELL * mat);
void matrix_ell_show(matrix_ELL * mat);

#endif
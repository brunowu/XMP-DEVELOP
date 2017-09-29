#include "mVecMatrix_complex.h"

void vector_divid(vector * v, complex * x){
	int length = vector_total(v);

	for(int i=0; i<length; i++){
		complex_divid((complex *)vector_get(v, i), x);
	}
}

void vector_multiple(vector * v, complex * x){
	int length = vector_total(v);

	for(int i=0; i<length; i++){
		complex_multiple((complex *)vector_get(v, i), x);
	}
}

void vector_add_vector(vector * v1, vector * v2, vector * newV){	
	if(vector_total(v1) != vector_total(v2)){
		printf("Not same dimension\n");
		return;
	}else{
		vector_init(newV, vector_total(v1));
		complex * items = (complex *)malloc(sizeof(complex) * vector_total(v1));
		for(int i=0; i<vector_total(v1); i++){
			complex_copy((items + i), (complex *)vector_get(v1, i));
			complex_add((items + i), (complex *)vector_get(v2, i));
			vector_add_duplicate(newV, (void *)(items + i));
		}
	}
}

void vector_reduce_vector(vector * v1, vector * v2, vector * newV){
	if(vector_total(v1) != vector_total(v2)){
		printf("Not same dimension\n");
		return;
	}else{
		vector_init(newV, vector_total(v1));
		complex * items = (complex *)malloc(sizeof(complex) * vector_total(v1));
		for(int i=0; i<vector_total(v1); i++){
			complex_copy((items + i), (complex *)vector_get(v1, i));
			complex_reduce((items + i), (complex *)vector_get(v2, i));
			vector_add_duplicate(newV, (void *)(items + i));
		}
	}
}

void vector_inner_produit(vector * v1, vector * v2, complex * h){
	if(vector_total(v1) != vector_total(v2)){
		printf("different dimension!\n");
	}else{
		int length = vector_total(v1);
		complex * item = (complex *)malloc(sizeof(complex));
		complex_init(h, 0, 0);
		for(int i=0; i<length; i++){
			complex_copy(item, (complex *)vector_get(v1, i));
			complex_add(item, (complex *)vector_get(v2, i));
			complex_add(h, item);
		}
	}
}

void vector_abs(vector * v, complex * h){
	vector_inner_produit(v, v, h);
	complex_sqrt(h);
}

void stock_vector_in_matrix(matrix * mat, vector * v, int ligne){
	if(matrix_get_colone(mat) > vector_total(v)){
		printf("Not matched colone!!!\n");
		return;
	}else{
		if(mat->ligne <= ligne){
			matrix_resize_ligne(mat, ligne + 1);
		}
		for(int i=0; i<vector_total(v); i++){
			if(i < matrix_total_ligne(mat, ligne)){
				matrix_set(mat, ligne, i, vector_get(v, i));
			}else{
				matrix_add_duplicate(mat, ligne, vector_get(v, i));
			}		
		}
	}
}

void stock_mvector_in_matrix(matrix * mat, void ** v, int ligne, int colone){
	if(matrix_get_colone(mat) > colone){
		printf("Not matched colone!!!\n");
		return;
	}else{
		if(mat->ligne <= ligne){
			if(mat->ligne == ligne){
				matrix_resize_ligne(mat, ligne + 1);
			}else{
				matrix_resize_ligne(mat, ligne);
			}	
		}

		for(int i=0; i<colone; i++){
			if(i < matrix_total_ligne(mat, ligne)){
				matrix_set(mat, ligne, i, v + i);
			}else{
				matrix_add_duplicate(mat, ligne, v + i);
			}
		}
	}
}

void matrix_get_vector(matrix * mat, vector * v, int ligne){
	vector_init(v, matrix_total_ligne(mat, ligne));
	for(int i=0; i<matrix_total_ligne(mat, ligne); i++){
		vector_add_duplicate(v, matrix_get(mat, ligne, i));
	}
}

void matrix_multiple_vector(matrix * mat, vector * v, vector * newV){
	int ligne = matrix_get_ligne(mat);
	int colone = matrix_get_colone(mat);

	if(colone != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		vector_init(newV, ligne);
		complex * items = (complex *)malloc(sizeof(complex) * ligne);
		complex * item = (complex *)malloc(sizeof(complex));

		for(int i=0; i<ligne; i++){
			complex_init(items + i, 0, 0);
			for(int j=0; j<colone; j++){
				complex_copy(item, (complex *)matrix_get(mat, i, j));
				complex_multiple(item, (complex *)vector_get(v, j));
				complex_add(items + i, item);
			}
		}
		for(int i=0; i<ligne; i++){
			vector_add_duplicate(newV, (void *)(items + i));
		}
		free(items); free(item);
	}	
}

//mat2 is a transpose of matrix
void matrix_multiple_matrix(matrix * mat1, matrix * mat2, matrix * newMat){
	vector v1, v2;
	matrix m;
	matrix_init(&m, matrix_get_ligne(mat2), matrix_get_ligne(mat1));
	for(int i=0; i<matrix_get_ligne(mat2); i++){
		matrix_get_vector(mat2, &v1, i);
		matrix_multiple_vector(mat1, &v1, &v2);
		stock_vector_in_matrix(&m, &v2, i);
		vector_free(&v1);
		vector_free(&v2);
	}
	matrix_transpose(&m, newMat);
	matrix_free(&m);
}

void upper_triangle_matrix_inverse(matrix * mat, matrix * mati){
	complex * item = (complex *)malloc(sizeof(complex));
	complex * index = (complex *)malloc(sizeof(complex));

	matrix_init(mati, matrix_get_ligne(mat), matrix_get_colone(mat));
	matrix_complete_ligne(mati);

	for(int k=0; k<matrix_get_ligne(mati); k++){
		if(((complex *)matrix_get(mat, k, k))->re == 0 && ((complex *)matrix_get(mat, k, k))->im == 0){
			printf("Matrix is not inversable!!!\n");
			exit(1);
		}else{
			complex_init(item, 1, 0);
			complex_divid(item, (complex *)matrix_get(mat, k, k));
			matrix_set(mati, k, k, (void *)item);
		}
		if(k > 0){
			for(int j=k-1; j>=0; j--){
				complex_init(item, 0, 0);
				for(int i=j+1; i<=k; i++){
					complex_copy(index, (complex *)matrix_get(mat, j, i));
					complex_multiple(index, (complex *)matrix_get(mati, i, k));
					complex_reduce(item, index);
				}
				if(((complex *)matrix_get(mat, j, j))->re == 0 && ((complex *)matrix_get(mat, j, j))->im == 0){
					printf("Matrix is not inversable!!!\n");
					exit(1);
				}else{
					complex_divid(item, (complex *)matrix_get(mat, j, j));
					matrix_set(mati, j, k, (void *)item);
				}
			}
		}
	}
	free(item); free(index);
}

/******************
Matrix Ellpack
*******************/
void matrix_ell_multiple_vector(matrix_ELL * mat, vector * v, vector * newV){
	if(mat->n_colone != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		vector_init(newV, vector_total(v));
		complex * index = (complex *)malloc(sizeof(complex));
		complex * items = (complex *)malloc(sizeof(complex) * mat->ligne);
		complex * item = (complex *)malloc(sizeof(complex));
		for(int i=0; i<mat->ligne; i++){
			complex_init(items + i, 0, 0);
			for(int j=0; j<mat->colone/2; j++){
				complex_copy(index, (complex *)matrix_ell_get(mat, i, j));
				if(index->re != -1 && index->im == 0){
					complex_copy(item, (complex *)matrix_ell_get(mat, i, mat->colone/2 + j));
					complex_multiple(item, (complex *)vector_get(v, index->re));
					complex_add(items + i, item);
				}	
			}
		}
		for(int i=0; i<mat->ligne; i++){
			vector_add_duplicate(newV, (void *)(items + i));
		}
		free(items); free(item); free(index);
	}	
}

//mat2 is a transpose of matrix
void matrix_ell_multiple_matrix(matrix_ELL * mat1, matrix * mat2, matrix * newMat){
	vector v1, v2;
	matrix m;
	matrix_init(&m, matrix_ell_get_ligne(mat1), matrix_get_ligne(mat2));
	for(int i=0; i<matrix_get_ligne(mat2); i++){
		matrix_get_vector(mat2, &v1, i);
		matrix_ell_multiple_vector(mat1, &v1, &v2);
		stock_vector_in_matrix(&m, &v2, i);
		vector_free(&v1);
		vector_free(&v2);
	}
	matrix_transpose(&m, newMat);
	matrix_free(&m);
}

void matrix_convert_matrix_ell(matrix * mat, matrix_ELL * matell){
	complex * item;
	int * count, * index;
	int k;
	count = (int *)malloc(sizeof(int));
	index = (int *)malloc(sizeof(int));
	item = (complex *)malloc(sizeof(complex));

	* count = 0;
	for(int i=0; i<mat->ligne; i++){
		* index = 0;
		for(int j=0; j<mat->colone; j++){
			if(!(((complex *)matrix_get(mat, i, j))->re == 0 && ((complex *)matrix_get(mat, i, j))->im == 0)){
				* index = * index + 1;
			}
		}
		//printf("max non-zero number of row %d is %d\n", i, * index);
		if(* index > * count){
			* count = * index;
		}
	}
	//printf("max non-zero number is %d\n", * count);
	matrix_ell_init(matell, mat->ligne, 2 * (* count), mat->colone);

	for(int i=0; i<mat->ligne; i++){
		k = 0;
		for(int j=0; j<mat->colone; j++){
			if(!(((complex *)matrix_get(mat, i, j))->re == 0 && ((complex *)matrix_get(mat, i, j))->im == 0)){
				complex_init(item, j, 0);
				matrix_ell_add_duplicate(matell, i, k, (void *)item);
				matrix_ell_add_duplicate(matell, i, (* count) + k, matrix_get(mat, i, j));
				k++;
			}
		}
	}
	free(count); free(index); free(item);
	matrix_ell_complete_ligne(matell);
}
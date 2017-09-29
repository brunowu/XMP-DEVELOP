#include "mVecMatrix.h"

void vector_divid(vector * v, double x){
	int length = vector_total(v);
	double * items;
	items = (double *)malloc(sizeof(double) * length);
	for(int i=0; i<length; i++){
		*(items + i) = (*(double *)vector_get(v, i))/x;
	}

	for(int i=0; i<length; i++){
		vector_set(v, i, (void *)(items + i));
	}
	free(items);
}

void vector_multiple(vector * v, double x){
	vector_divid(v, 1/x);
}

void vector_add_vector(vector * v1, vector * v2, vector * newV){
	
	if(vector_total(v1) != vector_total(v2)){
		printf("Not same dimension\n");
		return;
	}else{
		vector_init(newV, vector_total(v1));
		double * items = (double *)malloc(sizeof(double) * vector_total(v1));
		for(int i=0; i<vector_total(v1); i++){
			*(items + i) = *(double *)vector_get(v1, i) + *(double *)vector_get(v2, i);
			vector_add_duplicate(newV, (void *)(items + i));
		}
	}
}

void vector_reduce_vector(vector * v1, vector * v2, vector * newV){
	vector_divid(v2, -1);
	vector_add_vector(v1, v2, newV);
	vector_divid(v2, -1);
}

void vector_inner_produit(vector * v1, vector * v2, double * h){
	* h = 0;
	if(vector_total(v1) != vector_total(v2)){
		printf("different dimension!\n");
	}else{
		int length = vector_total(v1);
		* h = 0;
		for(int i=0; i<length; i++){
			* h += (*(double *)vector_get(v1, i)) * (*(double *)vector_get(v2, i));
		}
	}
}

void vector_abs(vector * v, double * h){
	vector_inner_produit(v, v, h);
	* h = sqrt(* h);
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
		double * items;
		items = (double *)malloc(sizeof(double) * ligne);
		for(int i=0; i<ligne; i++){
			*(items + i) = 0;
			for(int j=0; j<colone; j++){
				*(items + i) += (*(double *)matrix_get(mat, i, j)) * (*(double *)vector_get(v, j));
			}
		}
		for(int i=0; i<ligne; i++){
			vector_add_duplicate(newV, (void *)(items + i));
		}
		free(items);
		items = NULL;
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
	double * item;
	matrix_init(mati, matrix_get_ligne(mat), matrix_get_colone(mat));
	matrix_complete_ligne(mati);
	for(int k=0; k<matrix_get_ligne(mati); k++){
		item = (double *)malloc(sizeof(double));
		if(*(double *)matrix_get(mat, k, k) == 0){
			printf("Matrix is not inversable!!!\n");
			exit(1);
		}else{
			* item = 1 / (*(double *)matrix_get(mat, k, k));
			matrix_set(mati, k, k, (void *)item);
		}
		free(item);
		if(k > 0){
			for(int j=k-1; j>=0; j--){
				item = (double *)malloc(sizeof(double));
				* item = 0;
				for(int i=j+1; i<=k; i++){
					* item -= (*(double *)matrix_get(mat, j, i)) * (*(double *)matrix_get(mati, i, k)); 
				}
				if(*(double *)matrix_get(mat, j, j) == 0){
					printf("Matrix is not inversable!!!\n");
					exit(1);
				}else{
					* item = (* item) / (*(double *)matrix_get(mat, j, j));
					matrix_set(mati, j, k, (void *)item);
				}
				free(item);
			}
		}
	}
	item = NULL;
}

/******************
Matrix Ellpack
*******************/
void matrix_ell_multiple_vector(matrix_ELL * mat, vector * v, vector * newV){
	double index;
	if(mat->n_colone != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		vector_init(newV, vector_total(v));
		double * items;
		items = malloc(sizeof(double) * mat->ligne);
		for(int i=0; i<mat->ligne; i++){
			*(items + i) = 0;
			for(int j=0; j<mat->colone/2; j++){
				index = *(double *)matrix_ell_get(mat, i, j);
				if(index != (double)(-1)){
					*(items + i) += (*(double *)matrix_ell_get(mat, i, mat->colone/2 + j)) * (*(double *)vector_get(v, index));
				}	
			}
		}
		for(int i=0; i<mat->ligne; i++){
			vector_add_duplicate(newV, (void *)(items + i));
		}
		free(items);
		items = NULL;
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
	double * item;
	int * count, * index;
	int k;
	count = malloc(sizeof(int));
	index = malloc(sizeof(int));

	* count = 0;
	for(int i=0; i<mat->ligne; i++){
		* index = 0;
		for(int j=0; j<mat->colone; j++){
			if(*(double *)matrix_get(mat, i, j) != 0){
				* index = * index + 1;
			}
		}
		if(* index > * count){
			* count = * index;
		}
	}
	free(index); index = NULL;

	matrix_ell_init(matell, mat->ligne, 2 * (* count), mat->colone);

	for(int i=0; i<mat->ligne; i++){
		k = 0;
		for(int j=0; j<mat->colone; j++){
			if(*(double *)matrix_get(mat, i, j) != 0){
				item = malloc(sizeof(double));
				* item = (double)j;
				matrix_ell_add_duplicate(matell, i, k, (void *)item);
				matrix_ell_add_duplicate(matell, i, (* count) + k, matrix_get(mat, i, j));
				k++;
				free(item); item = NULL;
			}
		}
	}
	free(count); count = NULL;
	matrix_ell_complete_ligne(matell);
}






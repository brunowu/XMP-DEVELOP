#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constant_data.h"

//random_get function
double random_get(double min, double max);
int matrix_max_count(double ** m);

//main function
int main(void){
	srand(time(NULL));
	double ** mat;
	mat = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		mat[i] = malloc(sizeof(double) * COLS_NUM);
	}
	/**
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			if( j == i || j == i - 1){
				*(*(mat + i) + j) = random_get(-10, 10);
			}else if( j == i + 1){
				if(i > 0){
					*(*(mat + i) + j) = *(*(mat + i) + j - 2);
				}else if(i == 0){
					*(*(mat + i) + j) = random_get(-10, 10);
				}
			}else{
				*(*(mat + i) + j) = 0;
			}
		}
	}
	**/
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			if( j == i - 1){
				*(*(mat + i) + j) = -3;
			}else if( j == i ){
				*(*(mat + i) + j) = 5;
			}else if( j == i + 1){
				*(*(mat + i) + j) = 1;
			}else{
				*(*(mat + i) + j) = 0;
			}
		}
	}

	//write to mat_sample.txt
	FILE *f1;
	f1 = fopen("mat_10000_10000.txt", "wb+");
	int ret;
	for(int i=0; i<ROWS_NUM; i++){
		ret = fwrite(mat[i], sizeof(double), COLS_NUM, f1);
	}
	fclose(f1);
/**
	//matrix convert to matrix ellpack
	int count = 2 * matrix_max_count(mat);
	double ** mat_ell;
	mat_ell = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		mat_ell[i] = malloc(sizeof(double) * count);
	}
	int k;

	//initialize matrix ellpack
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<count; j++){
			*(*(mat_ell + i) + j) = -1;
		}
	}

	//copy data from matrix
	for(int i=0; i<ROWS_NUM; i++){
		k = 0;
		for(int j=0; j<COLS_NUM; j++){
			if(*(*(mat + i) + j) != 0){
				*(*(mat_ell + i) + k) = (double)j;
				*(*(mat_ell + i) + (count/2 + k)) = *(*(mat + i) + j);
				k++;
			}
		}
	}
	
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<count; j++){
			printf("%f ", *(*(mat_ell + i) + j));
		}
		printf("\n");
	}

	//write to mat_ell_sample.txt
	FILE *f2;
	f2 = fopen("mat_ell_gmres_simple_300_300.txt", "wb+");
	for(int i=0; i<ROWS_NUM; i++){
		ret = fwrite(mat_ell[i], sizeof(double), 2 * count, f2);
	}
	fclose(f2);
	free(mat_ell);
**/
	for(int i=0; i<ROWS_NUM; i++){
		free(mat[i]);
	}
	free(mat);
	printf("Well Done!!!\n");
}

double random_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

int matrix_max_count(double ** m){
	int * count, * index;
	count = malloc(sizeof(int));
	index = malloc(sizeof(int));

	* count = 0;
	for(int i=0; i<ROWS_NUM; i++){
		* index = 0;
		for(int j=0; j<COLS_NUM; j++){
			if(*(*(m + i) + j) != 0){
				* index = * index + 1;
			}
		}
		if(* index > * count){
			* count = * index;
		}
	}
	free(index); index = NULL;
	return (* count);
}
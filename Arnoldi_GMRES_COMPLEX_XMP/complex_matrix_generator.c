#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constant_data.h"
#include "complex_scalar.h"	

//main function
int main(void){
	srand(time(NULL));
	double ** mat;
	mat = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		mat[i] = malloc(sizeof(double) * 2 * COLS_NUM);
	}

	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<2 * COLS_NUM; j=j+2){
			if( j == 2*(i - 1)){
				*(*(mat + i) + j) = -3;
				*(*(mat + i) + j + 1) = 0;
			}else if( j == 2*i ){
				*(*(mat + i) + j) = 5;
				*(*(mat + i) + j + 1) = 0;
			}else if( j == 2*(i + 1)){
				*(*(mat + i) + j) = 1;
				*(*(mat + i) + j + 1) = 0;
			}else{
				*(*(mat + i) + j) = 0;
				*(*(mat + i) + j + 1) = 0;
			}
		}
	}

	//write to mat_sample.txt
	FILE *f1;
	f1 = fopen("mat_complex_10_10.txt", "wb+");
	int ret;
	for(int i=0; i<ROWS_NUM; i++){
		ret = fwrite(mat[i], sizeof(double), 2 * COLS_NUM, f1);
	}
	fclose(f1);

	for(int i=0; i<ROWS_NUM; i++){
		free(mat[i]);
	}
	free(mat);
	printf("Well Done!!!\n");
}
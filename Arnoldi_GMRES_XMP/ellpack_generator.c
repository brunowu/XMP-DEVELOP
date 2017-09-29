#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constant_data.h"

//random_get function
double random_double_get(double min, double max);
int random_int_get(double min, double max);

//generate the defined matrix of COLS_ELL_NUM = 6
void generate_defined_ellpack(double ** m);
//generate random matrix ellpack
void generate_random_ellpack(double ** m);

//main function
int main(void){
	srand(time(NULL));
	
	//allocate memory to matrix ellpack
	double ** mat_ell;
	mat_ell = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		mat_ell[i] = malloc(sizeof(double) * COLS_ELL_NUM);
	}

	generate_defined_ellpack(mat_ell);
	//generate_random_ellpack(mat_ell);

	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			printf("%f ", *(*(mat_ell + i) + j));
		}
		printf("\n");
	}

	//write to mat_ell_large_scale.txt
	FILE *f2;
	int ret;
	f2 = fopen("mat_ell_1000000_2002.txt", "wb+");
	for(int i=0; i<ROWS_NUM; i++){
		ret = fwrite(mat_ell[i], sizeof(double), COLS_ELL_NUM, f2);
	}
	fclose(f2);
}

void generate_defined_ellpack(double ** m){
	int index = (COLS_ELL_NUM/2 - 1)/2;
	//initialize matrix ellpack
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			*(*(m + i) + j) = -1;
		}
	}
	//generate data from matrix
	for(int i=0; i<ROWS_NUM; i++){
		if(i < index){
			for(int j=0; j<i; j++){
				*(*(m + i) + j) = j;
				*(*(m + i) + COLS_ELL_NUM/2 + j) = -3;
			}
			*(*(m + i) + i) = i;
			*(*(m + i) + COLS_ELL_NUM/2) = 5;
			for(int j=1; j<=index; j++){	
				*(*(m + i) + i + j) = i + j;
				*(*(m + i) + COLS_ELL_NUM/2 + i + j) = 1;	
			}		
		}else if(i > ROWS_NUM - 1 - index){
			for(int j=0; j<index; j++){
				*(*(m + i) + j) = i - index + j;
				*(*(m + i) + j + COLS_ELL_NUM/2) = -3;
			}	
			*(*(m + i) + index) = i;
			*(*(m + i) + index + COLS_ELL_NUM/2) = 5;
			for(int j=1; j<ROWS_NUM-i; j++){
				*(*(m + i) + index + j) = i + j;
				*(*(m + i) + index + j + COLS_ELL_NUM/2) = 1;
			}
		}else{
			for(int j=0; j<index; j++){
				*(*(m + i) + j) = i - index + j;
				*(*(m + i) + j + COLS_ELL_NUM/2) = -3;
			}
			*(*(m + i) + index) = i;
			*(*(m + i) + COLS_ELL_NUM/2 + index) = 5;
			for(int j=1; j<=index; j++){
				*(*(m + i) + index + j) = i + j;
				*(*(m + i) + COLS_ELL_NUM/2 + index + j) = 1;
			}	
		}
	}
}

void generate_random_ellpack(double ** m){
	int flag, index;
	for(int i=0; i<ROWS_NUM; i++){
		index = 0;
		flag = random_int_get(0, COLS_ELL_NUM/2);
		for(int j=0; j<=flag; j++){
			if(index < COLS_NUM - 1){
				index = random_int_get(index, COLS_NUM);
				*(*(m + i) + j) = index;
				*(*(m + i) + j + COLS_ELL_NUM/2) = random_double_get(-10, 10);
				index++;
			}else if(index == COLS_NUM - 1){
				*(*(m + i) + j) = index;
				*(*(m + i) + j + COLS_ELL_NUM/2) = random_double_get(-10, 10);
				index++;
			}else if(index > COLS_NUM - 1){
				*(*(m + i) + j) = -1;
				*(*(m + i) + j + COLS_ELL_NUM/2) = -1;
			}
		}
		for(int j=flag+1; j<COLS_ELL_NUM/2; j++){
			*(*(m + i) + j) = -1;
			*(*(m + i) + j + COLS_ELL_NUM/2) = -1;
		}
	}
}



double random_double_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

int random_int_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return (int)(min + (rand() / div));
}
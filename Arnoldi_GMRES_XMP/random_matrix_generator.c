#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constant_data.h"

//random_get function
double random_get(double min, double max);
void generate_random_matrix(double ** m);

int main(void){
	srand(time(NULL));
	double ** m;
	m = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m[i] = malloc(sizeof(double) * COLS_NUM);
	}

	generate_random_matrix(m);

	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			printf("%f ", *(*(m + i) + j));
		}
		printf("\n");
	}

	FILE *f1;
	int ret;
	f1 = fopen("mat_random_matrix_300_300.txt", "wb+");
	for(int i=0; i<ROWS_NUM; i++){
		ret = fwrite(m[i], sizeof(double), COLS_NUM, f1);
	}
	fclose(f1);
}

void generate_random_matrix(double ** m){
	//stock data in matrix
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			*(*(m + i) + j) = random_get(-10, 10);
		}
	}
}

double random_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return min + (rand() / div);
}
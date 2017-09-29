#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arnoldi_gmres.h"

#pragma xmp nodes p(NPES)
#pragma xmp template t(0:ROWS_NUM-1)
#pragma xmp distribute t(block) onto p

//matrix and vector for parallel computing
double mat[ROWS_NUM][COLS_NUM];
double mat_ell[ROWS_NUM][COLS_ELL_NUM];
double V[ROWS_NUM];
#pragma xmp align mat[i][*] with t(i)
#pragma xmp align mat_ell[ROWS_NUM][COLS_ELL_NUM] with t(i)
#pragma xmp align V[i] with t(i)

//global data
double m[ROWS_NUM][COLS_NUM];
double m_ell[ROWS_NUM][COLS_ELL_NUM];
//function which aide the routine
void readMatrix_matrix();
void readMatrix_ellpack();
void initialize_matrix(vector * v, matrix * matQ, matrix * matH);
void initialize_ellpack(vector * v, matrix * matQ, matrix * matH);
void Xmp_matrix_multiple_vector(vector * v);
void Xmp_ellpack_multiple_vector(vector * v);
void Xmp_vector_duplicate(double * v, vector * r);
void Xmp_matrix_arnoldi(vector * v, matrix * matQ, matrix * matH);
void Xmp_ellpack_arnoldi(vector * v, matrix * matQ, matrix * matH);

//main function of arnoldi algorithm
int main(void){
	double start_time, stop_time, elapsed_time;
	//Read matrix data
	readMatrix_matrix();
	//readMatrix_ellpack();
	#pragma xmp barrier

	#pragma xmp task on p(1)
{
	start_time = xmp_wtime();
}
	//Initialization of matrix and vector
	vector v;
	matrix matQ, matH, matT;
	initialize_matrix(&v, &matQ, &matH);
	//initialize_ellpack(&v, &matQ, &matH);
	#pragma xmp barrier

	//Arnoldi computing
	Xmp_matrix_arnoldi(&v, &matQ, &matH);
	//Xmp_ellpack_arnoldi(&v, &matQ, &matH);
	#pragma xmp barrier
	hessenberg_qr(&matH, &matT, IT_QR);
	#pragma xmp barrier

	#pragma xmp task on p(1)
{
	stop_time = xmp_wtime();
	elapsed_time = stop_time - start_time;

	printf("Total_time was %f seconds.\n", elapsed_time);
}
	vector_free(&v);
	matrix_free(&matQ);
	matrix_free(&matH);
	matrix_free(&matT);
}

/********************************
Arnoldi algorithm
*********************************/
void Xmp_matrix_arnoldi(vector * v, matrix * matQ, matrix * matH){
	vector r, r1;
	double * h = (double *)malloc(sizeof(double));
	double * v_p = (double *)malloc(sizeof(double) * ROWS_NUM);
	//calcul and stock of q1
	vector q1, q;
	vector_copy(v, &q1);
	vector_abs(v, h);
	vector_divid(&q1, * h);
	stock_vector_in_matrix(matQ, &q1, 0);
	free(h);
	vector_free(&q1);
	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_matrix_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &r);
		vector_free(&q);
		#pragma xmp barrier
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
			//vector_show(&q);
			//vector_show(&r);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			free(h);
			vector_free(&q); vector_free(&r1);
		}
		h = (double *)malloc(sizeof(double));
		vector_abs(&r, h);
		vector_divid(&r, * h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p);
			h = NULL; v_p = NULL;
			matrix_complete_ligne(matH);
			return;
			/**else{		
				matrix_free(matQ); matrix_init(matQ, restart_max, vector_total(v));
				matrix_free(matH); matrix_init(matH, restart_max, restart_max);				
				vector_free(v); free(h);
				arnoldi(mat, &r, matQ, matH, index_ligne, restart_max);
			}**/
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			vector_free(&r);
			free(h);
		}
	}
}

void Xmp_ellpack_arnoldi(vector * v, matrix * matQ, matrix * matH){
	vector r, r1;
	double * h = (double *)malloc(sizeof(double));
	double * v_p = (double *)malloc(sizeof(double) * ROWS_NUM);
	//calcul and stock of q1
	vector q1, q;
	vector_copy(v, &q1);
	vector_abs(v, h);
	vector_divid(&q1, * h);
	stock_vector_in_matrix(matQ, &q1, 0);
	free(h);
	vector_free(&q1);
	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_ellpack_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &r);
		vector_free(&q);
		#pragma xmp barrier
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
			//vector_show(&q);
			//vector_show(&r);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			free(h);
			vector_free(&q); vector_free(&r1);
		}
		h = (double *)malloc(sizeof(double));
		vector_abs(&r, h);
		vector_divid(&r, * h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p);
			h = NULL; v_p = NULL;
			matrix_complete_ligne(matH);
			return;
			/**else{		
				matrix_free(matQ); matrix_init(matQ, restart_max, vector_total(v));
				matrix_free(matH); matrix_init(matH, restart_max, restart_max);				
				vector_free(v); free(h);
				arnoldi(mat, &r, matQ, matH, index_ligne, restart_max);
			}**/
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			vector_free(&r);
			free(h);
		}
	}
}


/*****************************
Read Matrix
******************************/
void readMatrix_matrix()
{
	FILE * f1;
	double * temp = malloc(sizeof(double));

	f1 = fopen("mat_sample3.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			fread(temp, sizeof(double), 1, f1);
			m[i][j] = * temp;
		}
	}
}

void readMatrix_ellpack()
{
	FILE * f1;
	double * temp = malloc(sizeof(double));

	f1 = fopen("mat_ell_sample3.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			fread(temp, sizeof(double), 1, f1);
			m_ell[i][j] = * temp;
		}
	}
}

/*************************
Initialization of matrix
**************************/
void initialize_matrix(vector * v, matrix * matQ, matrix * matH){
	double a = 1;
	vector_init(v, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(v, (void *)&a);
	}
	matrix_init(matQ, restart_max, COLS_NUM);
	matrix_init(matH, restart_max, restart_max);

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			mat[i][j] = m[i][j];
		}
	}
}

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V[i] = 0;
	}
}
}

void initialize_ellpack(vector * v, matrix * matQ, matrix * matH){
	double a = 1;
	vector_init(v, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(v, (void *)&a);
	}
	matrix_init(matQ, restart_max, COLS_NUM);
	matrix_init(matH, restart_max, restart_max);

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			mat_ell[i][j] = m_ell[i][j];
		}
	}
}

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V[i] = 0;
	}
}
}

/*************************
Matrix Multiple Vector
**************************/
void Xmp_matrix_multiple_vector(vector * v){
	if(COLS_NUM != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		#pragma xmp loop on t(i)
{
		for(int i=0; i<ROWS_NUM; i++){
			V[i] = 0;
			for(int j=0; j<COLS_NUM; j++){
				V[i] += mat[i][j] * (*(double *)vector_get(v, j));
			}
		}
}	
	}
}

void Xmp_ellpack_multiple_vector(vector * v){
	if(COLS_NUM != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		#pragma xmp loop on t(i)
{
		for(int i=0; i<ROWS_NUM; i++){
			V[i] = 0;
			for(int j=0; j<COLS_ELL_NUM/2; j++){
				if(mat_ell[i][j] != -1){
					V[i] += mat_ell[i][COLS_ELL_NUM/2 + j] * (*(double *)vector_get(v, mat_ell[i][j]));
				}
			}
		}
}	
	}
}

void Xmp_vector_duplicate(double * v, vector * r){
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		*(v + i) = V[i];
	}
}
	vector_init(r, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(r, (void *)(v + i));
	}
}



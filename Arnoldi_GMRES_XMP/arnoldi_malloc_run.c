#include <time.h>
#include <xmp.h>
#include "arnoldi_gmres.h"
#include "constant_data.h"

#pragma xmp nodes p(NPES)
#pragma xmp template t(0:ROWS_NUM-1)
#pragma xmp distribute t(block) onto p

//matrix and vector for parallel computing
double (*mat)[COLS_NUM];
double (*mat_ell)[COLS_ELL_NUM];
double * V;
//double * vT;
//double * idx;

#pragma xmp align mat[i][*] with t(i)
#pragma xmp align mat_ell[i][*] with t(i)
#pragma xmp align V[i] with t(i)
//#pragma xmp align vT[i] with t(i)
//#pragma xmp align idx[i] with t(i)

//global data
double ** m;
double ** m_ell;
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
	//readMatrix_matrix();
	readMatrix_ellpack();
	#pragma xmp barrier

	#pragma xmp task on p(1)
{
	start_time = xmp_wtime();
}
	//Initialization of matrix and vector
	vector v;
	matrix matQ, matH, matT;
	//initialize_matrix(&v, &matQ, &matH);
	initialize_ellpack(&v, &matQ, &matH);
	#pragma xmp barrier

	//Arnoldi computing
	//Xmp_matrix_arnoldi(&v, &matQ, &matH);
	Xmp_ellpack_arnoldi(&v, &matQ, &matH);
	#pragma xmp barrier
	hessenberg_qr(&matH, &matT, IT_QR);
	#pragma xmp barrier

	#pragma xmp task on p(1)
{
	stop_time = xmp_wtime();
	elapsed_time = stop_time - start_time;
	matrix_show(&matT);
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
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&r, h);
		vector_divid(&r, * h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p);
			matrix_complete_ligne(matH);
			return;
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			vector_free(&r);
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
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&r, h);
		vector_divid(&r, * h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p);
			matrix_complete_ligne(matH);
			return;
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			vector_free(&r);
		}
	}
}


/*****************************
Read Matrix
******************************/
void readMatrix_matrix()
{
	//initialize matrix m
	m = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m[i] = malloc(sizeof(double) * COLS_NUM);
	}

	FILE * f1;
	double * temp = malloc(sizeof(double));

	f1 = fopen("mat_6000_6000.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			fread(temp, sizeof(double), 1, f1);
			*(*(m + i) + j) = * temp;
		}
	}
}

void readMatrix_ellpack()
{
	//Initialize of m_ell
	m_ell = malloc(sizeof(double *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m_ell[i] = malloc(sizeof(double) * COLS_ELL_NUM);
	}

	FILE * f1;
	double * temp = malloc(sizeof(double));

	f1 = fopen("mat_ell_100000_2002.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			fread(temp, sizeof(double), 1, f1);
			*(*(m_ell + i) + j) = * temp;
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

	mat = (double (*)[COLS_NUM])xmp_malloc(xmp_desc_of(mat), ROWS_NUM, COLS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			mat[i][j] = *(*(m + i) + j);
		}
	}
}

	V = (double *)xmp_malloc(xmp_desc_of(V), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V[i] = 0;
	}
}
/**
	vT = (double *)xmp_malloc(xmp_desc_of(vT), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		vT[i] = 0;
	}
}

	idx = (double *)xmp_malloc(xmp_desc_of(idx), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		idx[i] = i;
	}
}
**/
	for(int i=0; i<ROWS_NUM; i++){
		free(m[i]);
	}
	free(m);
}

void initialize_ellpack(vector * v, matrix * matQ, matrix * matH){
	double a = 1;
	vector_init(v, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(v, (void *)&a);
	}
	matrix_init(matQ, restart_max, COLS_NUM);
	matrix_init(matH, restart_max, restart_max);

	mat_ell = (double (*)[COLS_ELL_NUM])xmp_malloc(xmp_desc_of(mat_ell), ROWS_NUM, COLS_ELL_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			mat_ell[i][j] = *(*(m_ell + i) + j);
		}
	}
}
	V = (double *)xmp_malloc(xmp_desc_of(V), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V[i] = 0;
	}
}
/**
	vT = (double *)xmp_malloc(xmp_desc_of(vT), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		vT[i] = 0;
	}
}

	idx = (double *)xmp_malloc(xmp_desc_of(idx), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		idx[i] = i;
	}
}
**/
	for(int i=0; i<ROWS_NUM; i++){
		free(m_ell[i]);
	}
	free(m_ell);
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
/**	
	for(int i=0; i<ROWS_NUM; i++){
		*(v + i) = 0;
	}

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		*(v + i) = V[i];
	}
}

	for(int i=0; i<ROWS_NUM; i++){
		index = *(v + i);
		#pragma xmp barrier
		#pragma xmp reduction (+:index)
		*(v + i) = index;
	}
**/
	#pragma xmp gmove
{
	v[0:ROWS_NUM] = V[0:ROWS_NUM];
}
	vector_init(r, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(r, (void *)(v + i));
	}
}
/**
void Xmp_vector_gather(vector * r){
	xmp_gather(xmp_desc_of(vT),xmp_desc_of(V),xmp_desc_of(idx));

	#pragma xmp loop on t(i)
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		vT[i] = 0;
	}
}
}
**/
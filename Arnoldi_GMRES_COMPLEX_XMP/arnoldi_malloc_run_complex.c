#include <time.h>
#include <xmp.h>
#include "basic_operation_complex.h"
#include "constant_data.h"

#pragma xmp nodes p(NPES)
#pragma xmp template t(0:ROWS_NUM-1)
#pragma xmp distribute t(block) onto p

//matrix and vector for parallel computing
complex (*mat)[COLS_NUM];
complex (*mat_ell)[COLS_ELL_NUM];
double * V1;
double * V2;
//double * vT;
//double * idx;

#pragma xmp align mat[i][*] with t(i)
#pragma xmp align mat_ell[i][*] with t(i)
#pragma xmp align V1[i] with t(i)
#pragma xmp align V2[i] with t(i)
//#pragma xmp align vT[i] with t(i)
//#pragma xmp align idx[i] with t(i)

//global data
complex ** m;
complex ** m_ell;
//function which aide the routine
void readMatrix_matrix();
void readMatrix_ellpack();
void initialize_matrix(vector * v, matrix * matQ, matrix * matH);
void initialize_ellpack(vector * v, matrix * matQ, matrix * matH);
void Xmp_matrix_multiple_vector(vector * v);
void Xmp_ellpack_multiple_vector(vector * v);
void Xmp_vector_duplicate(double * v1, double * v2, vector * r);
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
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			complex_show(*(m + i) + j);
		}
		printf("\n");
	}
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
	//matrix_show(&matT);
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
	complex * h = (complex *)malloc(sizeof(complex));
	double * v_p1 = (double *)malloc(sizeof(double) * ROWS_NUM);
	double * v_p2 = (double *)malloc(sizeof(double) * ROWS_NUM);
	//calcul and stock of q1
	vector q1, q;
	vector_copy(v, &q1);
	vector_abs(v, h);
	vector_divid(&q1, h);
	stock_vector_in_matrix(matQ, &q1, 0);
	vector_free(&q1);
	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_matrix_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p1, v_p2, &r);
		vector_free(&q);
		#pragma xmp barrier
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&r, h);
		vector_divid(&r, h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p1); free(v_p2);
			matrix_complete_ligne(matH);
			return;
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			matrix_show(matQ);
			vector_free(&r);
		}
	}
}

void Xmp_ellpack_arnoldi(vector * v, matrix * matQ, matrix * matH){
	vector r, r1;
	complex * h = (complex *)malloc(sizeof(complex));
	double * v_p1 = (double *)malloc(sizeof(double) * ROWS_NUM);
	double * v_p2 = (double *)malloc(sizeof(double) * ROWS_NUM);
	//calcul and stock of q1
	vector q1, q;
	vector_copy(v, &q1);
	vector_abs(v, h);
	vector_divid(&q1, h);
	stock_vector_in_matrix(matQ, &q1, 0);
	vector_free(&q1);
	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_ellpack_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p1, v_p2, &r);
		vector_free(&q);
		#pragma xmp barrier
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, h);
			vector_reduce_vector(&r, &q, &r1);
			vector_free(&r);
			vector_duplicate(&r1, &r);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&r, h);
		vector_divid(&r, h);
		if(j == restart_max){
			vector_free(&r); vector_free(v);
			free(h); free(v_p1); free(v_p2);
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
	m = malloc(sizeof(complex *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m[i] = malloc(sizeof(complex) * COLS_NUM);
	}

	FILE * f1;
	complex * temp = malloc(sizeof(complex));
	double * flag = malloc(sizeof(double));

	f1 = fopen("mat_complex_10_10.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			fread(flag, sizeof(double), 1, f1); temp->re = * flag;
			fread(flag, sizeof(double), 1, f1); temp->im = * flag;
			complex_copy(*(m + i) + j, temp);
		}
	}
	free(temp); free(flag);
}

void readMatrix_ellpack()
{
	//Initialize of m_ell
	m_ell = malloc(sizeof(complex *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m_ell[i] = malloc(sizeof(complex) * COLS_ELL_NUM);
	}

	FILE * f1;
	complex * temp = malloc(sizeof(complex));
	double * flag = malloc(sizeof(double));

	f1 = fopen("mat_ell_100000_2002.txt", "rb");
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			fread(flag, sizeof(double), 1, f1); temp->re = * flag;
			fread(flag, sizeof(double), 1, f1); temp->im = * flag;
			complex_copy(*(m + i) + j, temp);
		}
	}
	free(temp); free(flag);
}

/*************************
Initialization of matrix
**************************/
void initialize_matrix(vector * v, matrix * matQ, matrix * matH){
	complex * a = (complex *)malloc(sizeof(complex));
	complex_init(a, 1, 0);
	vector_init(v, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(v, (void *)a);
	}
	matrix_init(matQ, restart_max, COLS_NUM);
	matrix_init(matH, restart_max, restart_max);

	mat = (complex (*)[COLS_NUM])xmp_malloc(xmp_desc_of(mat), ROWS_NUM, COLS_NUM);

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			complex_copy(*(mat + i) + j, *(m + i) + j);
		}
	}
}

	V1 = (double *)xmp_malloc(xmp_desc_of(V1), ROWS_NUM);
	V2 = (double *)xmp_malloc(xmp_desc_of(V2), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V1[i] = 0;
		V2[i] = 0;
	}
}

	for(int i=0; i<ROWS_NUM; i++){
		free(m[i]);
	}
	free(m); free(a);
}

void initialize_ellpack(vector * v, matrix * matQ, matrix * matH){
	complex * a = (complex *)malloc(sizeof(complex));
	complex_init(a, 1, 0);
	vector_init(v, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(v, (void *)&a);
	}
	matrix_init(matQ, restart_max, COLS_NUM);
	matrix_init(matH, restart_max, restart_max);

	mat_ell = (complex (*)[COLS_ELL_NUM])xmp_malloc(xmp_desc_of(mat_ell), ROWS_NUM, COLS_ELL_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			complex_copy(*(mat_ell + i) + j, *(m_ell + i) + j);
		}
	}
}
	V1 = (double *)xmp_malloc(xmp_desc_of(V1), ROWS_NUM);
	V2 = (double *)xmp_malloc(xmp_desc_of(V2), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V1[i] = 0;
		V2[i] = 0;
	}
}

	for(int i=0; i<ROWS_NUM; i++){
		free(m_ell[i]);
	}
	free(m_ell); free(a);
}

/*************************
Matrix Multiple Vector
**************************/
void Xmp_matrix_multiple_vector(vector * v){
	if(COLS_NUM != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		complex * a = (complex *)malloc(sizeof(complex));
		complex * flag = (complex *)malloc(sizeof(complex));
		complex_init(a, 0, 0);
		#pragma xmp loop on t(i)
{
		for(int i=0; i<ROWS_NUM; i++){
			V1[i] = 0; V2[i] = 0;
			for(int j=0; j<COLS_NUM; j++){
				complex_copy(flag, (complex *)vector_get(v, j));
				complex_multiple(flag, *(mat + i) + j);
				V1[i] += flag->re;
				V2[i] += flag->im;
			}
		}
}	
	free(a); free(flag);
	}
}

void Xmp_ellpack_multiple_vector(vector * v){
	if(COLS_NUM != vector_total(v)){
		printf("Not same dimension of matrix and vector");
		return;
	}else{
		complex * a = (complex *)malloc(sizeof(complex));
		complex * flag = (complex *)malloc(sizeof(complex));
		complex_init(a, 0, 0); 
		#pragma xmp loop on t(i)
{
		for(int i=0; i<ROWS_NUM; i++){
			V1[i] = 0; V2[i] = 0;
			for(int j=0; j<COLS_ELL_NUM/2; j++){
				if(!((*(mat_ell + i) + j)->re == -1 && (*(mat_ell + i) + j)->im == 0)){
					complex_copy(flag, (complex *)vector_get(v, (*(mat_ell + i) + j)->re));
					complex_multiple(flag, *(mat_ell + i) + COLS_ELL_NUM/2 + j);
					V1[i] += flag->re;
					V2[i] += flag->im;
				}
			}
		}
}	
	}
}

void Xmp_vector_duplicate(double * v1, double * v2, vector * r){
	#pragma xmp gmove
{
	v1[0:ROWS_NUM] = V1[0:ROWS_NUM];
}
	#pragma xmp gmove
{
	v2[0:ROWS_NUM] = V2[0:ROWS_NUM];
}
	vector_init(r, ROWS_NUM);
	complex * flag = (complex *)malloc(sizeof(complex));
	for(int i=0; i<ROWS_NUM; i++){
		complex_init(flag, v1[i], v2[i]);
		vector_add_duplicate(r, (void *)flag);
	}
}
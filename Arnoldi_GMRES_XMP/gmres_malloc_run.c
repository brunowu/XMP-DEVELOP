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

#pragma xmp align mat[i][*] with t(i)
#pragma xmp align mat_ell[i][*] with t(i)
#pragma xmp align V[i] with t(i)

//global data
double ** m;
double ** m_ell;
//function which aide the routine
void readMatrix_matrix();
void readMatrix_ellpack();
void initialize_matrix(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR);
void initialize_ellpack(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR);
void Xmp_matrix_multiple_vector(vector * v);
void Xmp_ellpack_multiple_vector(vector * v);
void Xmp_vector_duplicate(double * v, vector * r);
void Xmp_matrix_gmres(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int nb_tour);
void Xmp_ellpack_gmres(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int nb_tour);

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
	vector vx, vb;
	matrix matQ, matH, matT, omega, matR;
	initialize_matrix(&vx, &vb, &matQ, &matH, &omega, &matR);
	//initialize_ellpack(&vx, &vb, &matQ, &matH, &omega, &matR);
	#pragma xmp barrier

	//Gmres computing
	Xmp_matrix_gmres(&vx, &vb, &matQ, &matH, &omega, &matR, 0);
	//Xmp_ellpack_gmres(&vx, &vb, &matQ, &matH, &omega, &matR, 0);
	#pragma xmp barrier

	#pragma xmp task on p(1)
{
	stop_time = xmp_wtime();
	printf("vx capacity: %d, vb capacity: %d\n", vector_total(&vx), vector_total(&vb));
	elapsed_time = stop_time - start_time;
	printf("Total_time was %f seconds.\n", elapsed_time);
}
	vector_free(&vx);
	vector_free(&vb);
	matrix_free(&matQ);
	matrix_free(&matH);
	matrix_free(&matR);
	matrix_free(&omega);
}

/********************************
GMRES algorithm
*********************************/
void Xmp_matrix_gmres(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int nb_tour){
	#pragma xmp task on p(1)
{
	printf("Tour %d\n", nb_tour + 1);
}
	vector r, r1, rp, q, y, newx;
	matrix qT;
	double * h = (double *)malloc(sizeof(double));
	double * beta = (double *)malloc(sizeof(double));
	double * v_p = (double *)malloc(sizeof(double) * ROWS_NUM);

	//Calcul and stock q1
	#pragma xmp barrier
	Xmp_matrix_multiple_vector(vx);
	#pragma xmp barrier
	Xmp_vector_duplicate(v_p, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	* beta = * h;
	vector_divid(&r, * h);

	stock_vector_in_matrix(matQ, &r, 0);
	vector_free(&r);

	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_matrix_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &r);
		vector_duplicate(&r, &rp);
		vector_free(&q);
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&rp, &q, &r1);
			vector_free(&rp);
			vector_duplicate(&r1, &rp);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&rp, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&rp, * h);
		vector_free(&r);
		matrix_complete_ligne(matH);

		//linear_least_squares
		linear_least_squares(matH, omega, matR, &y, beta);

		matrix_transpose(matQ, &qT);
		//matrix_show(matQ); matrix_show(&qT); vector_show(&y);
		matrix_multiple_vector(&qT, &y, &q);
		matrix_free(&qT);
		vector_add_vector(vx, &q, &newx);
		vector_free(&q); vector_free(&y);
		#pragma xmp barrier
		Xmp_matrix_multiple_vector(&newx);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &q);
		vector_reduce_vector(vb, &q, &r1);
		vector_abs(&r1, h);
		vector_free(&q); vector_free(&r1); 
		#pragma xmp task on p(1)
{
		printf("residual is %.10f\n", * h);
}
		if( * h < 1e-10){
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}else if(* h > 1e10){
			printf("Divergence!!!\nQuit!!!\n");
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}
		if(j == restart_max){
			double * A = malloc(sizeof(double));
			* A = (double)1;
			nb_tour++;
			vector_free(&rp);
			free(h); free(beta); free(v_p);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
		    matrix_free(matQ); matrix_free(matH); matrix_free(omega); matrix_free(matR);
		    matrix_init(matQ, 1, 1); matrix_init(matH, 1, 1); matrix_init(matR, 1, 1); matrix_init(omega, 1, 1);
		    matrix_add_duplicate(omega, 0, (void *)A); 
		    free(A);
			Xmp_matrix_gmres(vx, vb, matQ, matH, omega, matR, nb_tour);
		}
		stock_vector_in_matrix(matQ, &rp, j);
		vector_free(&rp); vector_free(&newx);
	}
}

void Xmp_ellpack_gmres(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int nb_tour){
		#pragma xmp task on p(1)
{
	printf("Tour %d\n", nb_tour + 1);
}
	vector r, r1, rp, q, y, newx;
	matrix qT;
	double * h = (double *)malloc(sizeof(double));
	double * beta = (double *)malloc(sizeof(double));
	double * v_p = (double *)malloc(sizeof(double) * ROWS_NUM);

	//Calcul and stock q1
	#pragma xmp barrier
	Xmp_ellpack_multiple_vector(vx);
	#pragma xmp barrier
	Xmp_vector_duplicate(v_p, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	* beta = * h;

	vector_divid(&r, * h);
	stock_vector_in_matrix(matQ, &r, 0);
	vector_free(&r);

	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		#pragma xmp barrier
		Xmp_ellpack_multiple_vector(&q);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &r);
		vector_duplicate(&r, &rp);
		vector_free(&q);
		#pragma xmp barrier
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&rp, &q, &r1);
			vector_free(&rp);
			vector_duplicate(&r1, &rp);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&rp, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&rp, * h);
		vector_free(&r);
		matrix_complete_ligne(matH);

		//linear_least_squares
		linear_least_squares(matH, omega, matR, &y, beta);

		matrix_transpose(matQ, &qT);
		matrix_multiple_vector(&qT, &y, &q);
		//matrix_show(&qT); vector_show(&q);
		matrix_free(&qT);
		vector_add_vector(vx, &q, &newx);
		vector_free(&q); vector_free(&y);
		#pragma xmp barrier
		Xmp_ellpack_multiple_vector(&newx);
		#pragma xmp barrier
		Xmp_vector_duplicate(v_p, &q);
		vector_reduce_vector(vb, &q, &r1);
		vector_abs(&r1, h);
		vector_free(&q); vector_free(&r1); 
		printf("residual is %f\n", * h);
		if( * h < 3000){
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}
		if(j == restart_max){
			double * A = malloc(sizeof(double));
			* A = (double)1;
			nb_tour++;
			vector_free(&rp);
			free(h); free(beta); free(v_p);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
		    matrix_free(matQ); matrix_free(matH); matrix_free(omega); matrix_free(matR);
		    matrix_init(matQ, 1, 1); matrix_init(matH, 1, 1); matrix_init(matR, 1, 1); matrix_init(omega, 1, 1);
		    matrix_add_duplicate(omega, 0, (void *)A); 
		    free(A);
			Xmp_ellpack_gmres(vx, vb, matQ, matH, omega, matR, nb_tour);
		}
		stock_vector_in_matrix(matQ, &rp, j);
		vector_free(&rp); vector_free(&newx);
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

	f1 = fopen("mat_10000_10000.txt", "rb");
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
void initialize_matrix(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR){
	double a = 1;
	vector_init(vx, ROWS_NUM); vector_init(vb, ROWS_NUM);

	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(vx, (void *)&a);
		vector_add_duplicate(vb, (void *)&a);
	}
	
	matrix_init(matQ, 1, 1);
	matrix_init(matH, 1, 1);
	matrix_init(matR, 1, 1);
	matrix_init(omega, 1, 1);
	matrix_add_duplicate(omega, 0, (void *)&a);

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
	for(int i=0; i<ROWS_NUM; i++){
		free(m[i]);
	}
	free(m);
}

void initialize_ellpack(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR){
	double a = 1;
	vector_init(vx, ROWS_NUM); vector_init(vb, ROWS_NUM);
	
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(vx, (void *)&a);
		vector_add_duplicate(vb, (void *)&a);
	}

	matrix_init(matQ, 1, 1);
	matrix_init(matH, 1, 1);
	matrix_init(matR, 1, 1);
	matrix_init(omega, 1, 1);
	matrix_add_duplicate(omega, 0, (void *)&a);

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
	

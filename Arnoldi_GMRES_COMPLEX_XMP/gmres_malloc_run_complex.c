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
complex * V;

#pragma xmp align mat[i][*] with t(i)
#pragma xmp align mat_ell[i][*] with t(i)
#pragma xmp align V[i] with t(i)

//global data
complex ** m;
complex ** m_ell;
//function which aide the routine
void readMatrix_matrix();
void readMatrix_ellpack();
void initialize_matrix(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR);
void initialize_ellpack(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR);
void Xmp_matrix_multiple_vector(vector * v);
void Xmp_ellpack_multiple_vector(vector * v);
void Xmp_vector_duplicate(complex * v, vector * r);
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
	complex * h = (complex *)malloc(sizeof(complex));
	complex * beta = (complex *)malloc(sizeof(complex));
	complex * v_p = (complex *)malloc(sizeof(complex) * ROWS_NUM);

	//Calcul and stock q1
	#pragma xmp barrier
	Xmp_matrix_multiple_vector(vx);
	#pragma xmp barrier
	Xmp_vector_duplicate(v_p, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	complex_copy(beta, h);
	vector_divid(&r, h);

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
			vector_multiple(&q, h);
			vector_reduce_vector(&rp, &q, &r1);
			vector_free(&rp);
			vector_duplicate(&r1, &rp);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&rp, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&rp, h);
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
		printf("The residual is %f\n", complex_abs(h));
}
		if( complex_abs(h) < 1e-10){
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}else if(complex_abs(h) > 1e10){
			printf("Divergence!!!\nQuit!!!\n");
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}
		if(j == restart_max){
			complex * A = malloc(sizeof(complex));
			complex_init(A, 1, 0);
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
	complex * h = (complex *)malloc(sizeof(complex));
	complex * beta = (complex *)malloc(sizeof(complex));
	complex * v_p = (complex *)malloc(sizeof(complex) * ROWS_NUM);

	//Calcul and stock q1
	#pragma xmp barrier
	Xmp_ellpack_multiple_vector(vx);
	#pragma xmp barrier
	Xmp_vector_duplicate(v_p, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	complex_copy(beta, h);

	vector_divid(&r, h);
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
			vector_multiple(&q, h);
			vector_reduce_vector(&rp, &q, &r1);
			vector_free(&rp);
			vector_duplicate(&r1, &rp);
			vector_free(&q); vector_free(&r1);
		}
		vector_abs(&rp, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&rp, h);
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
		printf("residual is %f\n", complex_abs(h));
		if( complex_abs(h) < 3000){
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			return;
		}
		if(j == restart_max){
			complex * A = malloc(sizeof(complex));
			complex_init(A, 1, 0);
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
	m = malloc(sizeof(complex *) * ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		m[i] = malloc(sizeof(complex) * COLS_NUM);
	}

	FILE * f1;
	complex * temp = malloc(sizeof(complex));
	double * flag = malloc(sizeof(double));

	f1 = fopen("mat_complex_1000_1000.txt", "rb");
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
void initialize_matrix(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR){
	complex * a = (complex *)malloc(sizeof(complex));
	complex_init(a, 1, 0);
	vector_init(vx, ROWS_NUM); vector_init(vb, ROWS_NUM);

	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(vx, (void *)a);
		vector_add_duplicate(vb, (void *)a);
	}
	
	matrix_init(matQ, 1, 1);
	matrix_init(matH, 1, 1);
	matrix_init(matR, 1, 1);
	matrix_init(omega, 1, 1);
	matrix_add_duplicate(omega, 0, (void *)a);
	complex_init(a, 0, 0);

	mat = (complex (*)[COLS_NUM])xmp_malloc(xmp_desc_of(mat), ROWS_NUM, COLS_NUM);

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_NUM; j++){
			complex_copy(*(mat + i) + j, *(m + i) + j);
		}
	}
}

	V = (complex *)xmp_malloc(xmp_desc_of(V), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		complex_copy(V + i, a);
	}
}
	for(int i=0; i<ROWS_NUM; i++){
		free(m[i]);
	}
	free(m); free(a);
}

void initialize_ellpack(vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR){
	complex * a = (complex *)malloc(sizeof(complex));
	complex_init(a, 1, 0);
	vector_init(vx, ROWS_NUM); vector_init(vb, ROWS_NUM);
	
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(vx, (void *)a);
		vector_add_duplicate(vb, (void *)a);
	}

	matrix_init(matQ, 1, 1);
	matrix_init(matH, 1, 1);
	matrix_init(matR, 1, 1);
	matrix_init(omega, 1, 1);
	matrix_add_duplicate(omega, 0, (void *)&a);
	complex_init(a, 0, 0);

	mat_ell = (complex (*)[COLS_ELL_NUM])xmp_malloc(xmp_desc_of(mat_ell), ROWS_NUM, COLS_ELL_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		for(int j=0; j<COLS_ELL_NUM; j++){
			complex_copy(*(mat_ell + i) + j, *(m_ell + i) + j);
		}
	}
}

	V = (complex *)xmp_malloc(xmp_desc_of(V), ROWS_NUM);
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		complex_copy(V + i, a);
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
			complex_copy(V + i, a);
			for(int j=0; j<COLS_NUM; j++){
				complex_copy(flag, (complex *)vector_get(v, j));
				complex_multiple(flag, *(mat + i) + j);
				complex_add(V + i, flag);
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
			complex_copy(V + i, a);
			for(int j=0; j<COLS_ELL_NUM/2; j++){
				if(mat_ell[i][j] != -1){
					complex_copy(flag, (complex *)vector_get(v, mat_ell[i][j]));
					complex_multiple(flag, *(mat_ell + i) + COLS_ELL_NUM/2 + j);
					complex_add(V + i, flag);
				}
			}
		}
}	
	}
}

void Xmp_vector_duplicate(complex * v, vector * r){
	#pragma xmp gmove
{
	v[0:ROWS_NUM] = V[0:ROWS_NUM];
}
	vector_init(r, ROWS_NUM);
	for(int i=0; i<ROWS_NUM; i++){
		vector_add_duplicate(r, (void *)(v + i));
	}
}
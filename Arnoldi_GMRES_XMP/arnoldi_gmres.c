#include "arnoldi_gmres.h"

void arnoldi(matrix * mat, vector * v, matrix * matQ, matrix * matH){
	
	vector r, r1;
	double * h = (double *)malloc(sizeof(double));
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
		matrix_multiple_vector(mat, &q, &r);
		vector_free(&q);
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
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
			free(h);
			h = NULL;
			matrix_complete_ligne(matH);
			return;
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			stock_vector_in_matrix(matQ, &r, j);
			vector_free(&r);
			free(h); 
		}
	}
}

void gmres(matrix * mat, vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int index_ligne, double * beta){
	printf("Tour %d\n", index_ligne + 1);
	vector r, r1, rp, q, y, newx;
	matrix qT;
	double * h = (double *)malloc(sizeof(double));

	vector_init(&y, 1);
	matrix_multiple_vector(mat, vx, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	* beta = * h;

	vector_divid(&r, * h);
	stock_vector_in_matrix(matQ, &r, 0);
	free(h);
	vector_free(&r);

	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		matrix_multiple_vector(mat, &q, &r);
		vector_duplicate(&r, &rp);
		vector_free(&q); 
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&rp, &q, &r1);
			vector_free(&rp);
			vector_duplicate(&r1, &rp);
			free(h); 
			vector_free(&q); vector_free(&r1);
		}
		h = (double *)malloc(sizeof(double));
		vector_abs(&rp, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&rp, * h);
		free(h); vector_free(&r);
		matrix_complete_ligne(matH);

		//linear_least_squares
		linear_least_squares(matH, omega, matR, &y, beta);

		matrix_transpose(matQ, &qT);
		matrix_multiple_vector(&qT, &y, &q);
		//matrix_show(&qT); vector_show(&q);
		matrix_free(&qT);
		vector_add_vector(vx, &q, &newx);
		vector_free(&q); vector_free(&y);
		matrix_multiple_vector(mat, &newx, &q);
		vector_reduce_vector(vb, &q, &r1);
		h = (double *)malloc(sizeof(double));
		vector_abs(&r1, h);
		vector_free(&q); vector_free(&r1); 
		printf("residual is %f\n", * h);
		if( * h < 1e-1){
			vector_free(&rp);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			h = NULL; beta = NULL;
			return;
		}
		if(j == restart_max){
			double * A = malloc(sizeof(double));
			* A = (double)1;
			index_ligne++;
			vector_free(&rp);
			free(h); h = NULL;
			free(beta);
			beta = (double *)malloc(sizeof(double));
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
		    matrix_free(matQ); matrix_free(matH); matrix_free(omega); matrix_free(matR);
		    matrix_add_duplicate(omega, 0, (void *)A); 
		    free(A); A = NULL;
			gmres(mat, vx, vb, matQ, matH, omega, matR, index_ligne, beta);
		}
		stock_vector_in_matrix(matQ, &rp, j);
		vector_free(&rp); vector_free(&newx);
		free(h); 
	}
}

void linear_least_squares(matrix * H, matrix * omega, matrix * matR, vector * y, double * beta){
	double * c, * s, * ro, * sigma;
	vector v1, v2, vg;
	matrix newR, in_newR, newomega;
	//calcul for given rotation: c, s
	c = (double *)malloc(sizeof(double));
	s = (double *)malloc(sizeof(double));
	ro = (double *)malloc(sizeof(double));
	sigma = (double *)malloc(sizeof(double));

	matrix_get_vector(omega, &v1, matrix_get_ligne(omega) - 1);
	matrix_get_vector(H, &v2, matrix_get_ligne(H) - 1);
	vector_delete(&v2, vector_total(&v2) - 1);
	vector_inner_produit(&v1, &v2, ro);
	* sigma = *(double *)matrix_get(H, matrix_get_ligne(H) - 1, matrix_get_colone(H) - 1);
	givens(c, s, ro, sigma);
	vector_free(&v2);
	//printf("c: %f, s: %f, ro: %f, sigma: %f\n", * c, * s, * ro, * sigma);
	//calcul for new omega
	vector_multiple(&v1, * c);
	stock_vector_in_matrix(omega, &v1, matrix_get_ligne(omega) - 1);
	vector_divid(&v1, * c); vector_multiple(&v1, -(* s));
	stock_vector_in_matrix(omega, &v1, matrix_get_ligne(omega));
	matrix_add_duplicate(omega, matrix_get_ligne(omega) - 2, (void *)s);
	matrix_add_duplicate(omega, matrix_get_ligne(omega) - 1, (void *)c);
	matrix_complete_ligne(omega);
	vector_free(&v1);
	free(c); free(s); free(ro); free(sigma);
	c = NULL; s = NULL; ro = NULL; sigma = NULL;
	//calcul for new matR
	matrix_free(matR); 
	matrix_multiple_matrix(omega, H, matR);
	//matrix_show(omega); matrix_show(H); matrix_show(matR);
	//calcul for yn
	matrix_delete_ligne(matR, &newR, matR->ligne - 1);
	matrix_transpose(omega, &newomega);
	matrix_get_vector(&newomega, &vg, 0);
	vector_multiple(&vg, * beta);
	vector_delete(&vg, vector_total(&vg) - 1);
	upper_triangle_matrix_inverse(&newR, &in_newR);
	matrix_multiple_vector(&in_newR, &vg, y);
	matrix_free(&newR); matrix_free(&in_newR); matrix_free(&newomega);
	vector_free(&vg);
}

void hessenberg_qr(matrix * matH, matrix * matU, int num_iteration){
	double * c, * s, * m, * n;
	matrix_transpose(matH, matU);

	for(int k=1; k<=num_iteration; k++){
		for(int i=0; i<matrix_get_ligne(matU) - 1; i++){
			c = (double *)malloc(sizeof(double));
			s = (double *)malloc(sizeof(double));
			givens(c, s, (double *)matrix_get(matU, i, i), (double *)matrix_get(matU, i + 1, i));
			for(int j=i; j<matrix_get_ligne(matU); j++){		
				m = (double *)malloc(sizeof(double));
				n = (double *)malloc(sizeof(double));
				* m = ((* c) * (*(double *)matrix_get(matU, i, j))) + ((* s) * (*(double *)matrix_get(matU, i + 1, j)));
				* n = ((* c) * (*(double *)matrix_get(matU, i + 1, j))) - ((* s) * (*(double *)matrix_get(matU, i, j)));
				matrix_set(matU, i, j, (void *)m);
				matrix_set(matU, i + 1, j, (void *)n);
				free(m); free(n);
				m = NULL; n = NULL;
			}
			free(c); free(s);
			c = NULL; s = NULL;
		}
		for(int i=0; i<matrix_get_ligne(matU) - 1; i++){
			c = (double *)malloc(sizeof(double));
			s = (double *)malloc(sizeof(double));
			givens(c, s, (double *)matrix_get(matU, i, i), (double *)matrix_get(matU, i + 1, i));
			for(int j=0; j<=i + 1; j++){
				m = (double *)malloc(sizeof(double));
				n = (double *)malloc(sizeof(double));
				* m = (* c) * (*(double *)matrix_get(matU, j, i)) + (* s) * (*(double *)matrix_get(matU, j, i + 1));
				* n = (* c) * (*(double *)matrix_get(matU, j, i + 1)) - (* s) * (*(double *)matrix_get(matU, j, i));
				matrix_set(matU, j, i, (void *)m);
				matrix_set(matU, j, i + 1, (void *)n);
				free(m); free(n);
				m = NULL; n = NULL;
			}
			free(c); free(s);
			c = NULL; s = NULL;
		}
	}
}

static void givens(double * c, double * s, double * a, double * b){
	* c = (* a) / sqrt((* a) * (* a) + (* b) * (* b));
	* s = (* b) / sqrt((* a) * (* a) + (* b) * (* b));
}

double random_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

/***************
Matrix Ellpack
****************/
/**
void arnoldi_ell(matrix_ELL * mat, vector * v, matrix * matQ, matrix * matH, int index_ligne, int restart_max){	    
	printf("Tour %d\n", index_ligne + 1);
	vector r;
	double * h = (double *)malloc(sizeof(double));
	//calcul and stock of q1
	vector q1;
	vector_copy(v, &q1);
	vector_abs(v, h);
	vector_divid(&q1, * h);
	stock_vector_in_matrix(matQ, &q1, 0);
	free(h);
	vector_free(&q1);
	vector_free(v);

	for(int j=1; j<=restart_max; j++){
		vector q;
		matrix_get_vector(matQ, &q, j - 1);
		matrix_ell_multiple_vector(mat, &q, &r);
		vector_free(&q);
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector r1;
			vector_reduce_vector(&r, &q, &r1);
			vector_copy(&r1, &r);
			free(h);
			vector_free(&q);
			vector_free(&r1);
		}
		h = (double *)malloc(sizeof(double));
		vector_abs(&r, h);
		if(* h < 1){
			vector_free(&r);
			free(h);
			h = NULL;
			matrix_complete_ligne(matH);
			return;
		}else{
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_divid(&r, * h);
			free(h);
			h = NULL;
			stock_vector_in_matrix(matQ, &r, j);
			if(j == restart_max){
				index_ligne++;
				matrix_free(matQ); matrix_free(matH);
				arnoldi_ell(mat, &r, matQ, matH, index_ligne, restart_max);
			}
		}
	}
}

void gmres_ell(matrix_ELL * mat, vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int index_ligne, int restart_max, double * beta){
	printf("Tour %d\n", index_ligne + 1);
	vector r, r1, q, y, newx;
	matrix qT;
	double * h = (double *)malloc(sizeof(double));

	vector_init(&y);
	matrix_ell_multiple_vector(mat, vx, &q);
	vector_reduce_vector(vb, &q, &r);
	vector_free(&q);

	vector_abs(&r, h);
	* beta = * h;

	vector_divid(&r, * h);
	stock_vector_in_matrix(matQ, &r, 0);
	free(h);
	vector_free(&r);

	for(int j=1; j<=restart_max; j++){
		matrix_get_vector(matQ, &q, j - 1);
		matrix_ell_multiple_vector(mat, &q, &r);
		vector_free(&q);
		// Gram-Schmidt orthogonalization
		for(int i=1; i<=j; i++){
			matrix_get_vector(matQ, &q, i - 1);
			h = (double *)malloc(sizeof(double));
			vector_inner_produit(&q, &r, h);
			matrix_add_duplicate(matH, j - 1, (void *)h);
			vector_multiple(&q, * h);
			vector_reduce_vector(&r, &q, &r1);
			vector_duplicate(&r1, &r);
			free(h);
			vector_free(&q);
			vector_free(&r1);
		}
		h = (double *)malloc(sizeof(double));
		vector_abs(&r, h);
		matrix_add_duplicate(matH, j - 1, (void *)h);
		vector_divid(&r, * h);
		free(h);
		matrix_complete_ligne(matH);

		//linear_least_squares
		linear_least_squares(matH, omega, matR, &y, beta);

		matrix_transpose(matQ, &qT);
		matrix_multiple_vector(&qT, &y, &q);
		//matrix_show(&qT); vector_show(&q);
		matrix_free(&qT);
		vector_add_vector(vx, &q, &newx);
		vector_free(&q); vector_free(&y);
		matrix_ell_multiple_vector(mat, &newx, &q);
		vector_reduce_vector(vb, &q, &r1);
		h = (double *)malloc(sizeof(double));
		vector_abs(&r1, h);
		vector_free(&q); vector_free(&r1); 
		printf("residual is %f\n", * h);
		if( * h < 1e-5){
			vector_free(&r);
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
			free(h); free(beta);
			h = NULL; beta = NULL;
			return;
		}
		if(j == restart_max){
			double * A = malloc(sizeof(double));
			* A = (double)1;
			index_ligne++;
			vector_free(&r);
			free(h); h = NULL;
			free(beta);
			beta = (double *)malloc(sizeof(double));
			vector_free(vx);
		    vector_duplicate(&newx, vx);
		    vector_free(&newx);
		    matrix_free(matQ); matrix_free(matH); matrix_free(omega); matrix_free(matR);
		    matrix_add_duplicate(omega, 0, (void *)A); 
		    free(A); A = NULL;
			gmres_ell(mat, vx, vb, matQ, matH, omega, matR, index_ligne, restart_max, beta);
		}
		stock_vector_in_matrix(matQ, &r, j);
		vector_free(&r); vector_free(&newx);
		free(h); 
	}
}
**/



	

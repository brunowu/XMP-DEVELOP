#include "mVecMatrix.h"

#define ARNOLDI(mat, vec, matQ, matH) arnoldi(&mat, &v, &matQ, &matH)
//#define ARNOLDI_XMP(a, ligne, colone, vec, matQ, matH, restart_max) arnoldi_xmp(&&a, ligne, colone, &v, &matQ, &matH, restart_max)
#define GMRES(mat, vx, vb, matQ, matH, omega, matR, index_ligne, beta) gmres(&mat, &vx, &vb, &matQ, &matH, &omega, &matR, index_ligne, &beta)
#define HESSENBERG_QR(matH, matU, num_iteration) hessenberg_qr(&matH, &matU, num_iteration)
#define LINEAR_LEAST_SQUARES(H, omega, matR, y, beta) linear_least_squares(&omega, &matR, &y, &beta)
#define GIVENS(c, s, a, b) givens(&c, &s, &a, &b);
#define RANDOM_GET(min, max) random_get(min, max);

void arnoldi(matrix * mat, vector * v, matrix * matQ, matrix * matH);
//void arnoldi_xmp(double a[ROWS_NUM][COLS_NUM], vector * v, matrix * matQ, matrix * matH, int restart_max); 
void gmres(matrix * mat, vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int index_ligne, double * beta);
void linear_least_squares(matrix * H, matrix * omega, matrix * matR, vector * y, double * beta);
void hessenberg_qr(matrix * matH, matrix * matU, int num_iteration);
static void givens(double * c, double * s, double * a, double * b);
double random_get(double min, double max);

/*****************
Matrix Ellpack
******************/
//#define ARNOLDI_ELL(mat, vec, matQ, matH, index_ligne, restart_max) arnoldi_ELL(&mat, &v, &matQ, &matH, index_ligne, restart_max)
//#define GMRES_ELL(mat, vx, vb, matQ, matH, omega, matR, index_ligne, restart_max, beta) gmres(&mat, &vx, &vb, &matQ, &matH, &omega, &matR, index_ligne, restart_max, &beta)
//void arnoldi_ell(matrix_ELL * mat, vector * v, matrix * matQ, matrix * matH, int index_ligne, int restart_max);
//void gmres_ell(matrix_ELL * mat, vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int index_ligne, int restart_max, double * beta);

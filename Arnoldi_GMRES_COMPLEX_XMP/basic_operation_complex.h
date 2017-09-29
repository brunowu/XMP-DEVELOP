#include "mVecMatrix_complex.h"

#define ARNOLDI(mat, vec, matQ, matH) arnoldi(&mat, &v, &matQ, &matH)
#define GMRES(mat, vx, vb, matQ, matH, omega, matR, index_ligne, beta) gmres(&mat, &vx, &vb, &matQ, &matH, &omega, &matR, index_ligne, &beta)
#define HESSENBERG_QR(matH, matU, num_iteration) hessenberg_qr(&matH, &matU, num_iteration)
#define LINEAR_LEAST_SQUARES(H, omega, matR, y, beta) linear_least_squares(&omega, &matR, &y, &beta)
#define GIVENS(c, s, a, b) givens(&c, &s, &a, &b);
#define RANDOM_GET(min, max) random_get(min, max);

void arnoldi(matrix * mat, vector * v, matrix * matQ, matrix * matH);
void gmres(matrix * mat, vector * vx, vector * vb, matrix * matQ, matrix * matH, matrix * omega, matrix * matR, int index_ligne, complex * beta);
void linear_least_squares(matrix * H, matrix * omega, matrix * matR, vector * y, complex * beta);
void hessenberg_qr(matrix * matH, matrix * matU, int num_iteration);
static void givens(complex * c, complex * s, complex * a, complex * b);

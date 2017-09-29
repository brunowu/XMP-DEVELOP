#include <stdio.h>
#include <stdlib.h>
#include <xmp.h>
#include "arnoldi_gmres.h"

#pragma xmp nodes p(NPES)
#pragma xmp template t(0:ROWS_NUM-1)
#pragma xmp distribute t(block) onto p

double * V;
#pragma xmp align V[i] with t(i)

int main(void){
	V = (double *)xmp_malloc(xmp_desc_of(V), ROWS_NUM);
	double * A;
	double index;
	A = malloc(sizeof(double) * ROWS_NUM);

	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		V[i] = 1;
	}
}
	for(int i=0; i<ROWS_NUM; i++){
		A[i] = 0;
	}
	
	#pragma xmp loop on t(i)
{
	for(int i=0; i<ROWS_NUM; i++){
		A[i] = V[i];
	}
}

	for(int i=0; i<ROWS_NUM; i++){
		index = A[i];
		#pragma xmp barrier
		#pragma xmp reduction (+:index)
		A[i] = index;
	}

	#pragma xmp task on p(1)

{
	for(int i=0; i<ROWS_NUM; i++){
		printf("%f ", A[i]);
	}
}

}

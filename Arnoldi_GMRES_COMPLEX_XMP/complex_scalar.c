#include "complex_scalar.h"

/*************************************
Basic operations
*************************************/
void complex_init(complex * comp, double re, double im){
	comp->re = re;
	comp->im = im;
}

void complex_copy(complex * comp1, complex * comp2){
	comp1->re = comp2->re;
	comp1->im = comp2->im;
}

void complex_show(complex * comp){
	if(comp->im >=0){
		printf("%f+%fi ", comp->re, comp->im);
	}else{
		printf("%f%fi ", comp->re, comp->im);
	}
}
/*************************************
Arithmetic operations
**************************************/
void complex_add(complex * comp1, complex * comp2){
	comp1->re += comp2->re;
	comp1->im += comp2->im;
}

void complex_reduce(complex * comp1, complex * comp2){
	comp1->re -= comp2->re;
	comp1->im -= comp2->im;
}

void complex_multiple(complex * comp1, complex * comp2){
	double * re = malloc(sizeof(double));
	double * im = malloc(sizeof(double));
	* re = comp1->re * comp2->re - comp1->im * comp2->im;
	* im = comp1->re * comp2->im + comp1->im * comp2->re;
	comp1->re = * re;
	comp1->im = * im;
	free(re); free(im);
}

void complex_divid(complex * comp1, complex * comp2){
	double * re = malloc(sizeof(double));
	double * im = malloc(sizeof(double));
	* re = comp1->re * comp2->re + comp1->im * comp2->im;
	* im = comp1->im * comp2->re - comp1->re * comp2->im;
	comp1->re = * re;
	comp1->im = * im;
	free(re); free(im);

	double * reim = malloc(sizeof(double));
	* reim = comp2->re * comp2->re + comp2->im * comp2->im;
	comp1->re /= * reim;
	comp1->im /= * reim;

	free(reim);
}

void complex_sqrt(complex * comp){
	double * a = (double *)malloc(sizeof(double)); * a = comp->re;
	double * b = (double *)malloc(sizeof(double)); * b = comp->im;

	double * s = (double *)malloc(sizeof(double)); * s = sqrt((* a) * (* a) + (* b) * (* b));

	comp->re = sqrt(2)/2 * sqrt((* s) + * a);
	comp->im = sqrt(2)/2 * sqrt((* s) - * a);

	free(a); free(b); free(s);
}

double complex_abs(complex * comp){
	return sqrt(comp->re * comp->re + comp->im * comp->im);
}
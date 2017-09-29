#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef COMPLEX_H
#define COMPLEX_H

#define COMPLEX_INIT(comp, re, im) complex_init(&comp, re, im)
#define COMPLEX_COPY(comp1, comp2) complex_copy(&comp1, &comp2)
#define COMPLEX_SHOW(comp) complex_show(&comp)

#define COMPLEX_ADD(comp1, comp2) complex_add(&comp1, &comp2)
#define COMPLEX_REDUCE(comp1, comp2) complex_reduce(&comp1, &comp2)
#define COMPLEX_MULTIPLE(comp1, comp2) complex_multiple(&comp1, &comp2)
#define COMPLEX_DIVID(comp1, comp2) complex_divid(&comp1, &comp2)

#define COMPLEX_SQRT(comp) complex_sqrt(&comp)
#define COMPLEX_ABS(comp) complex_abs(&comp)

typedef struct complex {
	double re;
	double im;
} complex;

void complex_init(complex * comp, double re, double im);
void complex_copy(complex * comp1, complex * comp2);
void complex_show(complex * comp);

void complex_add(complex * comp1, complex * comp2);
void complex_reduce(complex * comp1, complex * comp2);
void complex_multiple(complex * comp1, complex * comp2);
void complex_divid(complex * comp1, complex * comp2);

void complex_sqrt(complex * comp);
double complex_abs(complex * comp);
#endif

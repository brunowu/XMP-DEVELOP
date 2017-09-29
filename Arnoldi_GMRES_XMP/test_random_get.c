#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int random_int_get(double min, double max);
int random_int_get(double min, double max){
	double range = max - min;
	double div = RAND_MAX / range;
	return (int)(min + (rand() / div));
}

int main(void){
	int k;
	for(int i=0; i<10000; i++){
		k = random_int_get(0, 2);
	//	if(k == 100 || k == 2){
			printf("%d ", k);
	//	}
	}
}
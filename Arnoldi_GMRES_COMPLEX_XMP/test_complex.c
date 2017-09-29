#include "basic_operation_complex.h"
#include "constant_data.h"

int main(void){
	complex c1, c2, c3, c4, c5;
	//complex * c = malloc(sizeof(complex));

	complex_init(&c1, -1, 1);
	c2 = c1;
	complex_copy(&c3, &c1);
	complex_show(&c2); complex_show(&c3);
	complex_add(&c2, &c1); complex_show(&c2);
	complex_reduce(&c2, &c1); complex_show(&c2);
	complex_multiple(&c3, &c1); complex_show(&c3);
	complex_divid(&c3, &c1); complex_show(&c3);
	complex_sqrt(&c3); complex_show(&c3);
	printf("%f\n", complex_abs(&c3));
}
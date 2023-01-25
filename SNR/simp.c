#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f (double x){
	return sqrt(x);
}

/* Composite Simpson's 3/8 rule accroding to wiki */
// Only use this formula if total number of subinternal is multiple of 3

// Evaluate the integral using function
double simp_f (double (*f)(double)) {

	// Properties of integration
	double lower, upper, step, val, ans;
	val = 0.0;
	lower = 0.0;
	upper = 10.0;
	int scale = 1;
	int size = 9 * scale; 			// Multiple of 3 ensured
	step = (upper - lower) / size;
	
	// generate the intergration interval	
	double x_axis[size + 1];
	int i;
	x_axis[0] = lower;
	for (i = 1; i < (size + 1); i++) {
		x_axis[i] = x_axis[i - 1] + step;
		//printf("i: %d, x_axis is %lf\n", i, x_axis[i]);
	}

	// Implementation of simpson's 3/8 rule
	for (i = 0; i < (size + 1); i++){
		if (i == 0 || i == size)
			val = val + (*f)(x_axis[i]);
		else if ((i % 3) != 0)
			val = val + 3.0 * (*f)(x_axis[i]);
		else if ((i % 3) == 0)
			val = val + 2.0 * (*f)(x_axis[i]);
	}

	ans = (3.0 * step / 8.0) * val;
	printf("ans is %.10lf\n", ans);
	return ans;
}

// Evaluate the integral by accepting pointer
double simp_xy(double *x, double *y, int size) {
	// Here size refer to total number of element in *x
	// Check if the array (size-1) is multiple of 3
	if (((size - 1) % 3) != 0){
		//printf("Array (size - 1) not multiple of 3! miniusing 1/2 data...\n");
		if (((size - 1) % 3) == 1)
			size = size - 1;
		else if (((size - 1) % 3) == 2)
			size = size - 2;
	}

	double val = 0.0;
	int i;
	// Implementation of simpson's 3/8 rule
	for (i = 0; i < size; i++){
		if (i == 0 || i == size - 1)
			val = val + *(y + i);
		else if ((i % 3) != 0)
			val = val + 3.0 * *(y + i);
		else if ((i % 3) == 0)
			val = val + 2.0 * *(y + i);
	}
	double step = (*(x + size - 1) - *x) / (size - 1);
	double ans;
	ans = (3.0 * step / 8.0) * val;
	//printf("ans is %.10lf\n", ans);
	return ans;	
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "setting.h"
#include "func_lib.h"

int main() {

	double *data;
	struct Pa p;
	char file[500] = "C0_D2_MA1.40_mA0.091_MB1.40_mB0.091_RA12.0_nA0.8_spin-2.0_e0.0_rp8.0_pA0.10_pB0.10.txt";
	char path[500] = "";
	int row;

	data = get_data(file, path, &row, 23);
	p = get_const(file);

	int j, i;
	double temp[23];
	double *acc;
	for (j = 0; j < row-1; j++) {

		for (i = 0; i < 22; i++) {
			temp[i] = *(data + j * 23 + i+1);
		}
		
		acc = eqn_orb(temp, p);
		printf("acc[0] is %e, acc[1] is %e\n", acc[0], acc[1]);

	}
}

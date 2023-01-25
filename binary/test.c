#include <stdio.h>
#include "func_lib.h"

struct Pa get_const(char *);

int main() {

	struct Pa p;
	char file[500] = "C0_D2_MA1.40_mA0.091_MB1.40_mB0.091_RA12.0_nA0.8_spin-2.0_e0.0_rp8.0_pA0.10_pB0.10.txt";

	p = get_const(file);
	printf("PNc: %d, PN1: %d, PN2: %d PN3: %d, PN25: %d\n",
			p.pn.PNc, p.pn.PN1, p.pn.PN2, p.pn.PN3, p.pn.PN25);

	printf("MA: %lf", p.M_A);
	return 0;
}
 

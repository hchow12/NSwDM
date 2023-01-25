#include <stdio.h>
#include "setting.h"

void evolve(struct Pa);

void evolve_BH(struct Pa);

/* m = mass of dark matter core 
 * M = mass of normal matter
 */

int main()
{
	struct Pa p; // Defined in setting.h header file

	p.M_A = 1.4; // Star A total mass (M+m) in (solar mass)
	p.p_A = 0.0; // Star A m/M ratio
	p.m_A = p.p_A / (1 + p.p_A); // Star A m/(M+m) ratio 

	p.M_B = 1.4; // Star B total mass (M+m) in (solar mass)
	p.p_B = 0.0; // Star B m/M ratio
	p.m_B = p.p_B / (1 + p.p_B); // Star B m/(M+m) ratio

	p.n_A = 0.8; // Star A polytropic index
	p.R_0_A = 12.0; // Star A radius without tide and dark matter core in unit of (km)
	p.f_R_A = -2.0;	// Irrotation star A as Riemann ellipsoid
	p.r = 8.0; // Initial orbital separation in unit of (R_0_A)

	/* Set to 1 to trigger the corresponding post-Newtonian terms, 0 otherwise.
	 * They modifythe post-Newtonian terms in evolution equation
	 * and initial orbital frequency
	 */
	p.pn.PNc = 1; // Conservation Post-Newtonian terms
	p.pn.PN1 = 1; // PN1 order
	p.pn.PN2 = 1; // PN2 order 
	p.pn.PN25 = 1; // Disspiative PN2.5 order
	p.pn.PN3 = 1; // PN3 order

	if (p.pn.PNc == 0 ) {
		p.pn.PN1 = 0;
		p.pn.PN2 = 0;
		p.pn.PN3 = 0;
	}

	p.dark.on = 0; // Set to 1 to trigger dark repulsion, 0 otherise (again affect the initial orbital frequency and evolution equation)
	p.dark.alpha = 0.01; // Dark repulsion coupling constant alpha prime
	p.dark.mV = 20.0; // Dark repulsion mediator range in km

	p.print_evolve = 'y'; // Set to 'n' to simiplify the program message, 'y' otherise.
	p.initial_only = 0; // Set to 1 to only solve for initial condition. 0 otherwise.

	if (p.pn.PN25 == 0) {
		printf("Warning: No dissipative PN terms, stopping the evolution...\n");
		p.initial_only = 1; // If no dissipative PN2.5 terms, produce only initial condition
	}

	// Uncomment and rewrite the following codes to loop through different parameters
	/*
	int len = 2;
	double n_list[len];
	int i;
	for (i=0; i < len; i++)
	{
		n_list[i] = 1.0 * (5 + i) / 10.0;
		p.n_A = n_list[i];
		printf("current setting: equal 1.4 solar mass, 0%%DM \t n = % .3lf\n", p.n_A);
		evolve(p);
	}
	*/

	 evolve(p); // evolve the binary and produce the data file, comment out this line if choose to loop through using above commented codes
}

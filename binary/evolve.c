#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_const_mksa.h>
#include "func_lib.h"
#define MAX_PT 4000000 // Set the maximum discrete point of solving the ODE be 4,000,000

void evolve(struct Pa pa)
{	
	time_t start = time(0); // Get the current real time
	time_t now; // Declare now variable of type time_t to get time later
	double c = GSL_CONST_MKSA_SPEED_OF_LIGHT; // Speed of light
	double mass_const = 1.3271244e20 / (c * c); // Convert unit of solar mass to geometric unit
	double *kapb;
	kapb = kap_b(pa.n_A); // Get polytropic constant kappa_n and b_n, as array of size 2
	struct Pa pb;
	pb = pa; // Copy the binary configuration setting to struct variable pb
	pb.n_B = pa.n_A;
	pb.c_n_A = 3.0 / (5.0 - pa.n_A); // Polytropic constant c_n
	pb.c_n_B = pb.c_n_A;
	pb.b_n_A = *(kapb + 1); // Polytropic constant b_n (due to DM core)
	pb.b_n_B = pb.b_n_A;
	pb.kap_A = *kapb; // Polytropic constant kappa_n (due to compressible ellipsoid)
	pb.kap_B = pb.kap_A;
	pb.f_R_B = pb.f_R_A;

	pb.M_A = pa.M_A * (1.0 - pa.m_A) * mass_const; // Star A normal matter mass (m)
	pb.m_A = pa.M_A * pa.m_A * mass_const; // Star A dark matter core mass (m)
	pb.M_B = pa.M_B * (1.0 - pa.m_B) * mass_const;
	pb.m_B = pa.M_B * pa.m_B * mass_const;
	pb.R_0_A = pb.R_0_A * 1000.0; // Star A reference radius R_0 (m)
	pb.R_0_B = pb.R_0_A * pow(pb.M_A / pb.M_B, (1.0 - pb.n_A) / (3.0 - pb.n_B)); // Star B reference radius (m)

	pb.dark.mV = 1.0 / (pa.dark.mV * 1000.0); // convert to unit (m^-1)

	struct Eqm eqm = initial(pb); // Obtain the equlibrium configuration for initial values
	// Store the initial value in array x[]
	double x[] = {0.0, eqm.a1_A
				   , 0.0
				   , eqm.a1_B
				   , 0.0
				   , eqm.a2_A
				   , 0.0
				   , eqm.a2_B
				   , 0.0
				   , eqm.a3_A
				   , 0.0
				   , eqm.a3_B
				   , eqm.orb
				   , 0.0
				   , eqm.orb
				   , 0
				   , -1.0 * pb.f_R_A * eqm.a1_A * eqm.a2_A * eqm.orb / (eqm.a1_A * eqm.a1_A + eqm.a2_A * eqm.a2_A)
				   , -1.0 * pb.f_R_B * eqm.a1_B * eqm.a2_B * eqm.orb / (eqm.a1_B * eqm.a1_B + eqm.a2_B * eqm.a2_B)
				   , 0.0
				   , eqm.r_phy
				   , eqm.orb
				   , 0.0};

	char *name;
	name = string_format(pa); // Obtain the name of the data file accroding to binary configuration

	FILE *f_write_log; // Declare a file pointer for writing "log.txt" file
	char log_file_name[500];
	sprintf(log_file_name, "HEAD_%s", name); // The log file is named as HEAD_data file name
	
	// Write the details of binary configuration into a file "log.txt"
	printf("Writing Configuration file...\n");
	f_write_log = fopen(log_file_name, "w+"); // Make the file pointer point to the log_file_name file
	fprintf(f_write_log, "%s\n", name);
	fprintf(f_write_log," Binary Configuration \n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Star A:\t\t\tStar B:\n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"Total Mass:\t        %.5lf\t        %.5lf (solar mass)\n", pa.M_A, pa.M_B);
	fprintf(f_write_log,"Normal Matter (NM) Mass:     %.5lf\t%.5lf (m)\n", pb.M_A, pb.M_B);
	fprintf(f_write_log,"Dark Matter (DM) Mass:       %.5lf\t%.5lf (m)\n", pb.m_A, pb.m_B);
	fprintf(f_write_log,"DM/Total Mass ratio:\t     %.5lf\t%.5lf \n", pb.m_A/(pb.M_A+pb.m_A), pb.m_B/(pb.M_B + pb.m_B));
	fprintf(f_write_log,"DM/NM Mass ratio:\t     %.5lf\t%.5lf \n", pb.m_A/pb.M_A, pb.m_B/pb.M_B);
	fprintf(f_write_log,"Binary Mass ratio:\t     %.5lf\n", pa.M_A/pa.M_B);
	fprintf(f_write_log,"Star's Reference radius: %.1lf (km)\n", pa.R_0_A);
	fprintf(f_write_log,"Star's Rotation number: %.1lf \n", pa.f_R_A);
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log," Orbital Configurations \n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"Conservative PN terms:\t\t%d\n", pa.pn.PNc);
	fprintf(f_write_log,"First Order PN:\t\t\t%d\n", pa.pn.PN1);
	fprintf(f_write_log,"Second Order PN:\t\t%d\n", pa.pn.PN2);
	fprintf(f_write_log,"Third Order PN:\t\t\t%d\n", pa.pn.PN3);
	fprintf(f_write_log,"Disspiative 2.5PN terms:\t%d\n", pa.pn.PN25);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Presence of Yukawa force:\t%d\n", pa.dark.on);
	fprintf(f_write_log,"Yukawa coupling number:\t\t%.5lf (repulsive)\n", pa.dark.alpha);
	fprintf(f_write_log,"Yukawa mediator mass range:\t%.5lf (km)\n", pa.dark.mV);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Initial Binary separation:\t%.5lf (km)\n", pa.r * pa.R_0_A);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log," Initial values \n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Star A:\t\t\tStar B:\n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"a1:    %.5lf\ta1:    %.5lf (m)\n", x[1], x[3]);
	fprintf(f_write_log,"a2:    %.5lf\ta2:    %.5lf (m)\n", x[5], x[7]);
	fprintf(f_write_log,"a3:    %.5lf\ta3:    %.5lf (m)\n", x[9], x[11]);
	fprintf(f_write_log,"spin:  %.5e\tspin:  %.5e (rad/m)\n", x[12], x[14]);
	fprintf(f_write_log,"vort:  %.5e\tvort:  %.5e (rad/m)\n", x[16], x[17]);
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log," Orbital parameters \n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"orbital separation:  %.5lf (m)\n", x[19]);
	fprintf(f_write_log,"orbital frequency:   %.5e (m)\n", x[20]);
	fprintf(f_write_log,"==========================END==========================\n");
	fclose(f_write_log);

	// Exit the program before solving the ODE
	if (pa.initial_only == 1) {
		printf("Exiting the program...\n");
		exit(0);
	}

	int dim;
	dim = sizeof(x) / sizeof(x[0]); // Obtain the dimension of the system of fist order ODE
	double t = 0; // set the initial time of evolution to 0 (m)
    double *t_array; // Declate the time array to store the time during evolution.
	double *y_array[dim]; // Declare a 2D array to store the variable a_i, r, omega, phase, etc. 
	t_array = (double *) malloc(MAX_PT * sizeof(double));
	double *p; // A 1-D array to store the first time derivative obtained by rk45 method
	double h = 100; // Minimum step size h = 100 (m)
	double eps = 1e-6; // Set the relative error between rk4 and rk5 to 1e-6, this parameter control the step size h
	int i = 0;
	int j;
	
	for (j = 0; j < dim; j++)
		y_array[j] = (double *) malloc(MAX_PT * sizeof(double));

	t_array[0] = 0; // set the initial time of evolution to 0 (m)
	for (j = 0; j < dim; j++)
		y_array[j][0] = x[j]; // Copy the initial value to the first row of 2D array
	
	while((x[1] + x[3]) < x[19]) // Stop the evolution when a_1 + a'_1 < r (orbital separation) 
	{	
		if (t_array[i]*1000/c > 1000) { // Prevent the evolution running longer than binary life time of 1000ms
			printf("Binary life exceed 1000ms! Stopping the evolution...\n");
			break;
		}
		
		if (h >= 100) // If the current step size is large enough
		{
			p = rk45(eqn, x, h, dim, pb); // Obtain the next point by rk45 function

			double err_array[dim]; // declare a 1D array to store the relative error (rk4 - rk5) of each equation in the ODE
			int l;
			for(l=0; l < dim; l++)
				err_array[l] = p[dim + l]; // Each relative error is stored in p[dim + l]
			while(max(err_array, dim) > eps) // If the maximum relative error is larger than the required error eps = 1e-6
			{	
				h = 0.9 * h * pow((eps / max(err_array, dim)), 1.0 / 5.0); // Reduce the step size by h_new = 0.9 * h_old * (eps/max)^(1/5)
				free(p); // Free the array for storing the next point and relative error
				p = rk45(eqn, x, h, dim, pb); // Redo the calculation using smaller step size
				for(l=0; l < dim; l++)
					err_array[l] = p[dim + l]; // Store the new relative error calculated from smaller step size, this end the while loop
			}

			for (j = 0; j < dim; j++)
				y_array[j][i + 1] = p[j]; // Store the next evolution point into the 2D y_array

			t = t + h; // Calculate the current time of evolution
			t_array[i + 1] = t; // Store the time value of the evolution point into t_array
			h = 0.9 * h * pow((eps / max(err_array, dim)), 1.0 / 5.0); // Enlarge the step size for next iteration
			free(p);
		}

		else if (h < 100) // If the current step size is too small
		{
			h = 100; // Set the step size to 100 (m) dispite the potential larger relative error
			p = rk45(eqn, x, h, dim, pb); // Obtain the next evolution point 
	
			for (j = 0; j < dim; j++)
				y_array[j][i + 1] = p[j]; // Store the next evolution point into y_array

			t = t + h; // Calculate the current time of evolution
			t_array[i + 1] = t; // Store the time value of evolution point into t_array
			
			free(p);
		}

		time(&now); // Get the current real time
		if ((i % 3000) == 0 && pa.print_evolve == 'y'){ // Print the evolution progress if print_evolve is set to 'y' in binary_conf.c
			printf("Real time = %ld (s)\tbinary time = %.3lf (ms)\tdistance = %.3lf\tstep = %.10lf\n",
				   now - start, t * 1000.0 / 3.0e8, y_array[19][i], h);
		}
		for (j = 0; j < dim; j++)
			x[j] = y_array[j][i + 1]; // Store the current evolution point into x[] to prepare for next iteration
		i++;
	}

	printf("\nSaving data to file...\n");

	FILE *fp;
	char naming[200]; // Store the directory+file_name for the output location of datafile
	if (pa.dark.on == 0)
		sprintf(naming, "data/dm/%s", name); // Store the data file into data/dm/ directory if no dark force
	else if (pa.dark.on == 1)
		sprintf(naming, "data/dark/%s", name); // Store the data file into data/dark/ directory if dark force is presence
	fp = fopen(naming, "w");
	int k = 0;
	// Writing the datafile from y_array and t_array
	while (k < i)
	{
		fprintf(fp, "%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",
		 t_array[k], y_array[0][k], y_array[1][k], y_array[2][k], y_array[3][k], y_array[4][k], y_array[5][k], y_array[6][k], y_array[7][k],
		 y_array[8][k], y_array[9][k], y_array[10][k], y_array[11][k], y_array[12][k], y_array[13][k], y_array[14][k], y_array[15][k],
		 y_array[16][k], y_array[17][k], y_array[18][k], y_array[19][k], y_array[20][k], y_array[21][k]);
		k++;
	}
	
	fclose(fp); // Close the data file

	free(t_array); // Free the t_array

	for (j = 0; j < dim; j++)
		free(y_array[j]); // Free the 2D y_array

	time(&now);
	printf("\nTotal time elapsed = %ld (s)\n", now - start);
	printf("\n*****************Program Done*****************\n");
}

/* This function evolve the BH-BH instead of NS-NS, 
 * the code here is largerly similar to the above NS-NS 
 * except the ODE dimension reduced to 4
 */
void evolve_BH(struct Pa pa)
{	
	time_t start = time(0);
	double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	double mass_const = 1.3271244e20 / (c * c);
	struct Pa pb;
	pb = pa;

	pb.M_A = pa.M_A * mass_const; 
	pb.M_B = pa.M_B * mass_const;
	pb.dark.mV = 1.0 / (pa.dark.mV * 1000.0);

	struct Eqm eqm = initial_BH(pb); // Obtain the initial value from intial_BH function
	time_t now;
	double x[] = {0.0, pb.r * pb.R_0_A * 1000, eqm.orb, 0.0};

	int dim;
	dim = sizeof(x) / sizeof(x[0]);

	char *name;
	name = string_format_BH(pa);

	FILE *f_write_log;
	char log_file_name[300];
	sprintf(log_file_name, "HEAD_%s", name);

	printf("Writing Configuration file...\n");
	f_write_log = fopen(log_file_name, "w+");	// Make the file pointer point to the log.txt file
	fprintf(f_write_log, "%s\n", name);
	fprintf(f_write_log," Binary Configuration \n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Star A:\t\t\tStar B:\n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"Individual Mass:\t        %.5lf\t        %.5lf (solar mass)\n", pa.M_A, pa.M_B);
	fprintf(f_write_log,"Binary Mass ratio:\t     %.5lf\n", pa.M_A/pa.M_B);
	fprintf(f_write_log,"Star's Reference radius: %.1lf (km)\n", pa.R_0_A);
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log," Orbital Configurations \n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"Conservative PN terms:\t\t%d\n", pa.pn.PNc);
	fprintf(f_write_log,"First Order PN:\t\t\t%d\n", pa.pn.PN1);
	fprintf(f_write_log,"Second Order PN:\t\t%d\n", pa.pn.PN2);
	fprintf(f_write_log,"Third Order PN:\t\t\t%d\n", pa.pn.PN3);
	fprintf(f_write_log,"Disspiative 2.5PN terms:\t%d\n", pa.pn.PN25);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Presence of Yukawa force:\t%d\n", pa.dark.on);
	fprintf(f_write_log,"Yukawa coupling number:\t\t%.5lf (repulsive)\n", pa.dark.alpha);
	fprintf(f_write_log,"Yukawa mediator mass range:\t%.5lf (km)\n", pa.dark.mV);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"Initial Binary separation:\t%.5lf (km)\n", pa.r * pa.R_0_A);
	fprintf(f_write_log,"\n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log," Initial values \n");
	fprintf(f_write_log,"==================================================\n");
	fprintf(f_write_log,"\n");
	fprintf(f_write_log," Orbital parameters \n");
	fprintf(f_write_log,"--------------------------------------------------\n");
	fprintf(f_write_log,"orbital separation:  %.5lf (m)\n", x[1]);
	fprintf(f_write_log,"orbital frequency:   %.5e (m)\n", x[2]);
	fprintf(f_write_log,"==========================END==========================\n");
	fclose(f_write_log);

	if (pa.initial_only == 1) {
		printf("Exiting the program...\n");
		exit(0);
	}

	double t = 0;
	double *t_array;
	double *y_array[dim];
	t_array = (double *) malloc(MAX_PT * sizeof(double));
	double *p;
	double h = 100;
	double eps = 1e-6;
	int i = 0;
	int j;
	
	for (j = 0; j < dim; j++)
		y_array[j] = (double *) malloc(MAX_PT * sizeof(double));

	t_array[0] = 0;

	for (j = 0; j < dim; j++)
		y_array[j][0] = x[j];
	
	while(x[1] > 2.5 * pb.R_0_A * 1000) // If the binary separation r > 2.5 R_0_A (defined in binary_conf.c)
	{	
		
		if (h >= 100)
		{
			p = rk45(eqn_BH, x, h, dim, pb);

			double err_array[dim];
			int l;
			for(l=0; l < dim; l++)
				err_array[l] = p[dim + l];
			while(max(err_array, dim) > eps)
			{	
				h = 0.9 * h * pow((eps / max(err_array, dim)), 1.0 / 5.0);
				free(p);
				p = rk45(eqn_BH, x, h, dim, pb);
				for(l=0; l < dim; l++)
					err_array[l] = p[dim + l];
			}

			for (j = 0; j < dim; j++)
				y_array[j][i + 1] = p[j];

			t = t + h;
			t_array[i + 1] = t;
			h = 0.9 * h * pow((eps / max(err_array, dim)), 1.0 / 5.0);
			free(p);
		}

		else if (h < 100)
		{
			h = 100;
			p = rk45(eqn_BH, x, h, dim, pb);

			for (j = 0; j < dim; j++)
				y_array[j][i + 1] = p[j];

			t = t + h;
			t_array[i + 1] = t;
			
			free(p);
		}

		time(&now);
		if ((i % 3000) == 0 && pa.print_evolve == 'y'){
			printf("Real time = %ld (s)\tbinary time = %.3lf (ms)\tdistance = %.3lf\tstep = %.10lf\n",
				   now - start, t * 1000.0 / 3.0e8, y_array[1][i], h);
		}
		for (j = 0; j < dim; j++)
			x[j] = y_array[j][i + 1];
		i++;
		
	}

	printf("\nSaving data to file...\n");

	FILE * fp;
	char naming[200];
	if (pa.dark.on == 0)
		sprintf(naming, "data/bh/no_dark/%s", name);
	else if (pa.dark.on == 1)
		sprintf(naming, "data/bh/dark/%s", name);
	fp = fopen(naming, "w");
	int k = 0;
	while (k < i)
	{
		fprintf(fp, "%.15lf %.15lf %.15lf %.15lf %.15lf\n", t_array[k], y_array[0][k], y_array[1][k], y_array[2][k], y_array[3][k]);
		k ++;
	}
	
	fclose(fp);

	free(t_array);

	for (j = 0; j < dim; j++)
		free(y_array[j]);

	time(&now);
	printf("\nTotal time elapsed = %ld (s)\n", now - start);
	printf("\n*****************Program Done*****************\n");
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "func_lib.h"
#define MAX_ROW 2000000 // Set the maximum row to read from data file be 2,000,000

int main(int argc, char *argv[]) {

	double *data;
	struct Pa p;
	double *acc;
	double temp[23];
	double *t, *ra, *oa, *r, *rv, *o, *ov;
	double *hpc, *hc, *hp; 

	t = malloc(MAX_ROW * sizeof(double)); // time (m)
	ra = malloc(MAX_ROW * sizeof(double)); // second time derivative of r 
	oa = malloc(MAX_ROW * sizeof(double)); // second time derivative of theta (orbital phase)
	r = malloc(MAX_ROW * sizeof(double)); // orbital separation r
	rv = malloc(MAX_ROW * sizeof(double)); // first time derivative of r
	o = malloc(MAX_ROW * sizeof(double)); // orbital phase theta
	ov = malloc(MAX_ROW * sizeof(double)); // first time derivative of theta
	hc = malloc(MAX_ROW * sizeof(double)); // h_x GW polarization mode
	hp = malloc(MAX_ROW * sizeof(double)); // h_+ GW polarization mode

	// Getting raw evolution data
	char file[500];
	char path[500] = "";
	int row;
	strcpy(file, argv[1]);
	printf("Producing GW file for %s....\n", file);

	data = get_data(file, path, &row, 23); // Reading data from datafile and store it in data
	if (data == NULL) {
		printf("Error encounter in reading file, exiting GW....");
		exit(1);
	}
	p = get_const(file); // Getting binary parameters from file name

	int i, j;

	for (j = 0; j < row - 1; j++)
	{
		for (i = 0; i < 22; i++) {
			t[j] = *(data + j * 23); // Extracting time data from datafile
			temp[i] = *(data + j * 23 + i+1); // Extracting orbital evolution data from datafile
		}

		acc = eqn_orb(temp, p); // Use the data from data file calculation orbital acceleration d^2r/dt^2 and d^2 theta/dt^2
		//printf("acc[0] is %e, acc[1] is %e\n", acc[0], acc[1]);
		r[j] = temp[19];
		rv[j] = temp[18];
		ra[j] = acc[0];
		o[j] = temp[21];
		ov[j] = temp[20];
		oa[j] = acc[1];

		hpc = TT_gauge(r[j], rv[j], ra[j], o[j], ov[j], oa[j]); // Calculation of GW polarization mode

		hp[j] = *hpc; // h_+ polarization mode data
		hc[j] = *(hpc + 1); // h_x polarization mode data

		free(hpc);
		free(acc);		
	}

	// Store the GW data into "Output.txt" file
	FILE *fp;
	char gw_file[400] = "Output.txt";
	fp = fopen(gw_file, "w");
	for (i = 0; i < row - 1; i++)
		fprintf(fp, "%.15lf\t%.15lf\t%.15lf\n", t[i], hp[i], hc[i]);
	fclose(fp);
	printf("THE GW file saved as %s\n", gw_file);
	
	free(t);
	free(ra);
	free(oa);
	free(r);
	free(rv);
	free(o);
	free(ov);
	free(hc);
	free(hp);
	free(data);

	return 0;
}

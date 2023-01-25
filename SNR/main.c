#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_const_mksa.h>

int line_counter(char*, char*);

double *get_data(char*, char*, int*, int);

double simp_xy(double*, double*, int);

int main(int argc, char *argv[]) {

	// Input Checking
	if (argc < 2) {
		printf(" This program need F_GW file as argument!\n");
		exit(1);
	}

	// Obtain the detector Noise PSD data from aLIGODesign.txt
	double *noise;
	char file[30] = "aLIGODesign.txt";
	char path[30] = "";
	unsigned int i, len;
	noise = get_data(file, path, &len, 2);

	// Obtain the fourier transformed GW data from F_GW file
	double *f_gw;
	char gw_file[500]; 
	strcpy(gw_file, argv[1]);
	char gw_path[500] = "";
	unsigned int gw_len;
	f_gw = get_data(gw_file, gw_path, &gw_len, 5);

	// Extraction of x, y noise PSD data
	double *x = (double *) malloc(len * sizeof(double));
	double *y = (double *) malloc(len * sizeof(double));
	for (i = 0; i < len; i++) {
		*(x + i) = *(noise + 2 * i);
		*(y + i) = *(noise + 2 * i + 1);
	}
	free(noise);

	// Extraction of f1, f2, A1, A2 data from F_GW data file (frequency and amplitude)
	double *f1 = (double *) malloc(gw_len * sizeof(double));
	double *f2 = (double *) malloc(gw_len * sizeof(double));
	double *A1 = (double *) malloc(gw_len * sizeof(double));
	double *A2 = (double *) malloc(gw_len * sizeof(double));
	for (i = 0; i < gw_len; i++) {
		*(f1 + i) = *(f_gw + 5 * i + 1);
		*(A1 + i) = *(f_gw + 5 * i + 2);
		*(f2 + i) = *(f_gw + 5 * i + 3);
		*(A2 + i) = *(f_gw + 5 * i + 4);
	}
	free(f_gw);

	// Interpolation function from gsl
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, len);
	gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, gw_len);
	gsl_spline *spline3 = gsl_spline_alloc(gsl_interp_cspline, gw_len);
	
	gsl_spline_init(spline, x, y, len);
	gsl_spline_init(spline2, f1, A1, gw_len);
	gsl_spline_init(spline3, f2, A2, gw_len);

	free(x);
	free(y);
	free(f1);
	free(f2);
	free(A1);
	free(A2);

	// Initialzing the interpolation x_axis
	double min, max, step;
	min = 5.0; // Minimum frequency
	max = 4900.0; // Maxiumum frequency
	step = 1.0; // frequency step
	size_t N = (size_t) ((max - min) / step);

	// Printing the interpolated noise psd and prepare the x/y array to integrate 
	double xi; // Frequency
	double yi; // Detector Noise PSD
	double a1i; // Amplitude For Fourier transformed GW h_+/x
	double a2i; // Amplitude For Fourier transformed GW h_+/x
	double x_integrate[N]; // Variable of integration (frequnecy) for SNR calculation
	double y_integrate[N]; // Value of integrand for SNR calculation
	double signal_psd; // GW signal PSD data

	// File saving variable declaration for saving the GW signal PSD
	char *datafile, buffer[50];
	char psd_path[400] = "signal_psd/";
	int fd, size;
	datafile = (char *) malloc(400);
	strcpy(datafile, "PSD_");
	strcat(datafile, gw_file);
	strcat(psd_path, datafile);

	// Opening file to save PSD data
	fd = open(psd_path, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
	free(datafile);

	// Writing GW signal PSD data and computing the SNR integrand
	i = 0;
	for (xi = min; xi < max; xi += step) {
		yi = gsl_spline_eval(spline, xi, acc);
		a1i = gsl_spline_eval(spline2, xi, acc);
		a2i = gsl_spline_eval(spline3, xi, acc);
		//printf("%g %g %g %g\n", xi, yi, a1i, a2i);
		x_integrate[i] = xi;
		//printf("x_integrate[%d] = %g\n", i, x_integrate[i]);
		y_integrate[i] = (a1i * a1i + a2i * a2i) / (yi * yi);
		
		signal_psd = sqrt(4 * (a1i * a1i + a2i * a2i) * xi); // This value need to be mutiplied by 2 * reduced mass / (distance to source) 
		
		size = sprintf(buffer, "%.10lf %.10lf\n", xi, signal_psd); // Writing GW signal PSD data
		write(fd, buffer, size);
		i++;
	}

	close(fd);

	// Evaluation of SNR
	double snr;
	snr = simp_xy(x_integrate, y_integrate, N);
	printf("%e\n", snr);


	gsl_spline_free(spline);
	gsl_spline_free(spline2);
	gsl_spline_free(spline3);
	gsl_interp_accel_free(acc);

	return 0;
}


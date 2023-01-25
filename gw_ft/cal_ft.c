#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_const_mksa.h>
#include "func_lib.h"

double signal(double);

int main(int arg, char *args[]) {
	
	// Reading raw gw data file (h+ or h_\times)
	double *data;
	double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;

	// Read file name from command line argument ./a.out <file_name>
	char path[500] = "";
	char file[500];
	strcpy(file, args[1]);
	unsigned int i, j, len; // len is the number of row of the data file
	data = get_data(file, path, &len, 3);

	// Extracting each column from data
	double x[len], y1[len], y2[len];
	for (i = 0; i < len; i++) {
		x[i] = *(data + 3 * i) / c; // time data, converted to second
		y1[i] = *(data + 3 * i + 1); // h_+ data, dimensionless
		y2[i] = *(data + 3 * i + 2); // h_x data, dimensionless
	}	
	free(data);

	// Using the interpolation function from gsl
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, len);
	gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, len);

	gsl_spline_init(spline, x, y1, len); // Interpolation of the t-h_+ graph
	gsl_spline_init(spline2, x, y2, len); // Interpolation of the t-h_x graph
	

	// T = sampling period: interval between sampling points
	// output array is complex, hence 2*N
	// Bigger N = bigger resolution in frequency (more digit in Hz)
	// Defalt N = 10000
	int N = 10000;
	double T, in[N], out[2 * N], mag[N];
	
	// Initializing input array with sampling period T = 0.00031 second
	// With N = 1000, this is equivalent to signal duration of 0.31s
	// smaller T = larger maximum output frequency
	T = x[len - 1] / (double) N; 
	

	// To do DFT on two array and store them in external array
	int col;
	int count = N/2 + 1;
	double amp1[count], freq1[count], amp2[count], freq2[count];
	for (col = 0; col < 2; col++) { // col = 0: Fourier Transform on h_+ array
		if (col == 0) {	            // col = 1: Fourier Transform on h_x array 
			for (i = 0; i < N; i++) 
				in[i] = gsl_spline_eval (spline, T * i, acc);
		}
		else if (col == 1) {	
			for (i = 0; i < N; i++) 
				in[i] = gsl_spline_eval (spline2, T * i, acc);
		}

		// Deployment of Forward DFT (Discrete Fourier Transform) formula
		double re = 0;
		double im = 0;
		int k;
		for (i = 0; i < N; i++) {
			for (k = 0; k < N; k++) {
				re = re + in[k] * cos(2.0 * M_PI * i * k / (float) N); // Real part of fourier transformed h_+/x 
				im = im + in[k] * sin(2.0 * M_PI * i * k / (float) N) * -1.0; // Imaginary part of fourier transformed h_+/x
			}
			out[2 * i] = re; // Store the real part in out
			out[2 * i + 1] = im; // Also store the Imaginary part in out
			mag[i] = sqrt(re * re + im * im); // Magnitude of fourier transformed h_+/x

			// Set re and im back to 0 to prepare for next iteration of i
			re = im = 0;
		}

		/* Printing calcuated value before and after
		for (i = 0; i < N; i++) {
			printf("index: %d input: %e DFT (Magnitude): %e \n", i, in[i], mag[i]); 
		}*/

		//printf("============================================\n");

		// Extraction of the correct value of magnitude
		double amp[N];
		for (i = 0; i < N; i++) {
			if (i == 0 || i == (N / 2))
				amp[i] = mag[i] / N;
			else
				amp[i] = mag[i] / (N / 2);
		}

		// Extraction of the frequency in Hz
		double freq[N/2 + 1];
		for (i = 0; i <= N/2; i++) {
			freq[i] =  i / (N * T);
		}

		// Storing amp and freq in another array (freq12), and free the gsl function
		if (col == 0) { 
			for (i = 0; i <= N/2; i++) {
				amp1[i] = amp[i]; // amp1 store the amplitude of fourier transformed h_+
				freq1[i] = freq[i]; // freq1 store the frequency axis of fourier transformed h_+
			}
			gsl_spline_free(spline);
		}
		else if (col == 1) { 
			for (i = 0; i <= N/2; i++) {
				amp2[i] = amp[i]; // amp2 store the amplitude of fourier transformed h_x
				freq2[i] = freq[i]; // freq2 store the frequency axis of fourier transformed h_x
			}
			gsl_spline_free(spline2);
		}
	}	
	gsl_interp_accel_free(acc);

	// Saving Amplitude and frequency to file
	FILE *fp;
	char ft_file[100] = "FT_Output.txt";
	fp = fopen(ft_file, "w");
	for (i = 0; i <= N/2; i ++) {
		fprintf(fp, "%d %g %g %g %g\n", i, freq1[i], amp1[i], freq2[i], amp2[i]);
	}
	fclose(fp);
	return 0;
}

/* signal function is only for testing the Fourier transform */
double signal(double t) {
	
	double y1, y2, y;
	y1 = 5.0 + 10.0 * cos(2.0 * M_PI * 0.1 * t - M_PI / 2.0);
	y2 = 3.0 * cos(2.0 * M_PI * 5.0 * t);
	y = y1 + y2;
	
	/*double y;
	if (t <= 10 && t > 0)
		y = 1.0;
	else 
		y = 0.0;*/
	return y;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Return the total line number of data file
int line_counter(char file[], char path[]) {

	FILE *fp;
	int count = 0;
	char c;

	strcat(path, file);
	fp = fopen(path, "r");

	for (c = getc(fp); c != EOF; c = getc(fp))
		if (c == '\n')
			count = count + 1;

	fclose(fp);

	return count - 1;
}

// Return the data stored in data file as pointer
double *get_data (char file[], char path[], int* lenn, int col) {

	double *p;
	unsigned int i, j;
	int len;
	len = line_counter(file, path);
	*lenn = len;
	
	// Initialize the 2d array (actually 1D array)
	p = (double *) malloc(len * col * sizeof(double));

	FILE *fp = fopen(path, "r");
	// Store the file content in 2d array
	for (i = 0; i < len; i++) {
		for (j = 0; j < col; j++) {
			fscanf(fp, "%lf", p + i * col + j);
		}
	}
	fclose(fp);

	return p;  
}

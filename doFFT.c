#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#define MIN(a, b) (((a)<(b))?(a):(b))

void bstr2darr(char *bstr, double darr[], int *len, int n) {
	
	int i = 0;
	
	while ((bstr[i] != '\0') && (i < n)) {
		/* printf("Bit[%d] = %c, ", i, bstr[i]); */
		if (bstr[i] == '0')
			darr[i] = 0;
		else if (bstr[i] == '1')
			darr[i] = 1;
		/* printf("Number[%d] = %1.0lf, ", i, darr[i]); */
		i++;
	}
	/*printf("\n");*/
	
	*len = i;
	
	while (i < n) {
		
		darr[i] = 0;
		i++;
	}
}

void darr2bstr(char **bstr, double darr[], int len) {
	
	int i;
	
	*bstr = (char *) malloc(sizeof(char)*(len+1));
	
	for (i = 0; i < len; i++)
		if (darr[i] > 0)
			(*bstr)[i] = '1';
		else 
			(*bstr)[i] = '0';
	(*bstr)[len] = '\0';
	
}

void print_darr(double darr[], int n) {

	int i;
	printf("[ %1.0lf", darr[0]);
	for (i=1; i<n; i++)
		printf(", %1.0lf", darr[i]);
	printf("]\r\n");
}

double min_amplitude(fftw_complex cfs[], int n) {

	int i = 0;
	int k = (n >> 1) + 1;
	
	double m = cfs[0][0]*cfs[0][0] +cfs[0][1]*cfs[0][1];
	
	for (i=1; i < k; i++)
		m = MIN(m, cfs[i][0]*cfs[i][0] + cfs[i][1]*cfs[i][1]);
		
	return sqrt(m);

}


int main(int argc, char *argv[]) {

	char *bstr;
	int n, l;
	
	double *in, min_amp;
	fftw_complex *out;
	fftw_plan dft_plan;
	
	sscanf(argv[1], "%d", &n);
	bstr = argv[2];
	
	in = (double *) fftw_malloc(sizeof(double)*n);
	out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*n);
	
	dft_plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
	
	bstr2darr(bstr, in, &l, n);
	
	fftw_execute(dft_plan);
	min_amp = min_amplitude(out, n);
	
	printf("DFT Minimal amplitude: %1.10lf\r\n", min_amp);
	
	fftw_destroy_plan(dft_plan);
	fftw_free(in);
	fftw_free(out);
	
	/* printf("Converted to array of length %d: ", l); */
	/* print_darr(darr, n); */
	/* darr2bstr(&bstr2, darr, l); */
	/*printf("Converted back to string: %s\r\n", bstr2); */
}

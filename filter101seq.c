#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

void init_FFT(int len, int N, double **in, fftw_complex **out, fftw_plan *plan) {
	
	int i;
	
	*in = (double *) fftw_malloc(N * sizeof(double));
	*out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
	
	*plan = fftw_plan_dft_r2c_1d(N, *in, *out, FFTW_MEASURE);
	
	(*in)[0] = 1;
	for (i=1; i < N-1; i++)
		(*in)[i] = 0;
	(*in)[len-1] = 1;
}

void clear_FFT(double *in, fftw_complex *out, fftw_plan plan) {

	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);

}

int find_first_set(int n) {
	int pos = 0;
	while ((n>0) && ((n%2) == 0)) {
		n >>= 1;
		pos++;
	}
	return pos;	
}

void print_seq(int len, double seq[]) {
	
	int i;
	for (i=0; i<len; i++)
		printf("%1.0lf", seq[i]);
	
}

void do_filter(int len, int K, double sq_min, long *total, long *count, double in[], fftw_complex out[], fftw_plan plan) {
	
	int i;
	double val2;
	
	fftw_execute(plan);
	
	double min2 = out[0][0] * out[0][0] + out[0][1]*out[0][1];
	
	i = 1;
	
	while ((i < K) && (min2 >= sq_min)) {
	
		val2 = out[i][0] * out[i][0] + out[i][1]*out[i][1];
		if (val2 < min2)
			min2 = val2;
		i++;
	};
	
	if (min2 >= sq_min) {
		
		print_seq(len, in);
		printf("   %1.3lf\r\n", sqrt(min2));
		(*count)++;
	}
	(*total)++;
}

int main(int argc, char *argv[])
{

	int len, k, n, shift, half, last, i, j, pos_i, pos_j;
	
	long tot, cnt;
	
	double t, minsq;
	
	clock_t start, end;
	
	double *seq;
	fftw_complex *dft;
	fftw_plan pln;
	
	if (argc != 4)
	{
		printf("Wrong no. of command line args!\r\n");
		return 0;
	}
	
	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &len);
	sscanf(argv[3], "%lf", &t);
	
	printf("#Samples N=%d, length=%d, threshold t=%1.3f\n", n, len, t);
	printf("#Sequences and their FFT minima:\n");
	
	if (len < 4)
	{
		printf("#Length must be >= 4, got length=%d instead!\r\n", len);
		return 0;
	
	} else if (len > n)
	{
		printf("#n must be >= lenght=%d, got n=%d instead!\r\n", len, n);
		return 0;
	} else if (len > 64)
	{
		printf("#l=%d is too large: must be <=64!\r\n", n);
		return 0;
	}
	
	start = clock();
	
	init_FFT(len, n, &seq, &dft, &pln);
	
	//printf("#Initialized input array: ");
	//print_seq(len, seq);
	//printf("\n");
	
	k = (n>>1)+1;
	
	half = (len-2) >> 1;
	shift = (len % 2)+1;
	last = (1 << half)-1;
	
	tot = 0;
	cnt = 0;
	
	minsq = t*t;
	
	i = 0; j = 1;
	
	if (len % 2 == 0)
		
		while (i < last) {
	
			while (j <= last) {
		
				pos_j = find_first_set(j);
				seq[half+pos_j+1] = 1 - seq[half+pos_j+1];
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				j++;
		
				}
			
			i++;
			if (i == last)
				break;
				
			pos_i = find_first_set(i);
			seq[half-pos_i]=1-seq[half-pos_i];
		
			while (j > i+2) {
		
				j--;
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				seq[half+pos_j+1] = 1 - seq[half+pos_j+1];
				pos_j = find_first_set(j);
		
			}
		
			do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
			i++;
			pos_i = find_first_set(i);
			seq[half-pos_i]=1-seq[half-pos_i];
		
		}
	else
		while (i < last) {
	
			while (j <= last) {
		
				pos_j = find_first_set(j);
				seq[half+pos_j+2] = 1 - seq[half+pos_j+2];
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				seq[half+1] = 1 - seq[half+1];
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				j++;
		
			}
			
			i++;
			if (i == last)
				break;
				
			pos_i = find_first_set(i);
			seq[half-pos_i]=1-seq[half-pos_i];
		
			while (j > i+2) {
		
				j--;
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				seq[half+1] = 1 - seq[half+1];
				do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
				seq[half+pos_j+2] = 1 - seq[half+pos_j+2];
				pos_j = find_first_set(j);
		
			}
		
			do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
			seq[half+1] = 1 - seq[half+1];
			do_filter(len, k, minsq, &tot, &cnt, seq, dft, pln);
	
			i++;
			pos_i = find_first_set(i);
			seq[half-pos_i]=1-seq[half-pos_i];
		
		}
	
	printf("#Total: %ld\r\n", tot);
	printf("#Left:  %ld\r\n", cnt);
	
	clear_FFT(seq, dft, pln);
	
	end = clock();
	printf("#Time:  %lf sec.\r\n", (double) (end - start)/ CLOCKS_PER_SEC);
	
	return 0;
}
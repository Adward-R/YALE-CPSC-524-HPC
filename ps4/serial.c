#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>

#define MAX_N_BODY 50000
#define TIME_INTERVAL 128
#define G 1.0

void readdata(char filename[], int *N, int *K, double *DT, double *mass, double *x, double *v) {
	FILE *fp;
	int i;

	if ((fp = fopen(filename, "r")) == NULL) {
		printf("Data cannot be read!\n");
		exit(1);
	}
	fscanf(fp, "%d", N);
	fscanf(fp, "%d", K);
	fscanf(fp, "%lf", DT);
	for (i = 0; i < *N; ++ i) { fscanf(fp, "%lf", mass+i); }
	for (i = 0; i < *N; ++ i) { fscanf(fp, "%lf %lf %lf", x+i*3, x+i*3+1, x+i*3+2); }
	for (i = 0; i < *N; ++ i) { fscanf(fp, "%lf %lf %lf", v+i*3, v+i*3+1, v+i*3+2); }
	fclose(fp);
}

void report(int t, int N, double DT, double mass[], double x[], double v[]) {
	double cx[3], cv[3];
	double msum = .0;
	int i, k;

	for (k = 0; k < 3; ++ k) cx[k] = cv[k] = .0;
	for (i = 0; i < N; ++ i) {
		msum += mass[i];
		for (k = 0; k < 3; ++ k) {
			cx[k] += mass[i] * x[i*3 + k];
			cv[k] += v[i*3 + k];
		}
	}
	for (k = 0; k < 3; ++ k) {
		// printf("cx[%d] = %lf\n", k, cx[k]);
		cx[k] /= msum;
		cv[k] /= N;
	}

	printf("\nConditions after timestep %d (time = %lf) :\n\n", t, t * DT);
	printf("\tCenter of Mass:\t(%lf, %lf, %lf)\n", cx[0], cx[1], cx[2]);
	printf("\tAverage Velocity:\t(%lf, %lf, %lf)\n", cv[0], cv[1], cv[2]);
}

int main(int argc, char **argv) {
	int N;  // number of bodies
	int K;  // number of time steps
	double DT;  // time step size
	double mass[MAX_N_BODY];
	double x[MAX_N_BODY][3];  // position in each dimension
	double v[MAX_N_BODY][3];  // velocity in each dimension
	double a[MAX_N_BODY][3];  // acceleration in each dimension
	double F[3];  // acceleration in each dimension

	int i, j, k, t;
	double tmp, r2;

	readdata("data/testdata1", &N, &K, &DT, mass, x, v);
	// printf("%d %d %lf\n", N, K, DT);
	// for (i = 0; i < N; ++ i) printf("%lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);

	for (t = 0; t < K; ++ t) {
		if (t % TIME_INTERVAL == 0) report(t, N, DT, mass, x, v);
		for (i = 0; i < N; ++ i) {  // for each body
			F[0] = F[1] = F[2] = .0;
			for (j = 0; j < N; ++ j) {
				if (j == i) continue;
				r2 = pow(x[i][0] - x[j][0], 2) + pow(x[i][1] - x[j][1], 2) + pow(x[i][2] - x[j][2], 2);
				if (r2 > 25.0) continue;
				tmp = (G * mass[i] * mass[j]) / (r2 * sqrt(r2));  // quantative variable
				for (k = 0; k < 3; ++ k) F[k] += tmp * (x[j][k] - x[i][k]);  // vector variable
			}
			for (k = 0; k < 3; ++ k) a[i][k] = F[k] / mass[k];
			for (k = 0; k < 3; ++ k) v[i][k] += a[i][k] * (DT / 2.0);
			for (k = 0; k < 3; ++ k) x[i][k] += v[i][k] * DT;
		}
	}
	report(K, N, DT, mass, x, v);

	return 0;
}
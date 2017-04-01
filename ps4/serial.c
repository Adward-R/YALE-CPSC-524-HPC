#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>

#define MAX_N_BODY 50000
#define TIME_INTERVAL 128
#define G 1.0
#define N_DIM 3
#define CUTOFF_SQR 25.0
#define EPSILON 0.0001

typedef struct body_t {
	double mass;
	double x[N_DIM];
	double v[N_DIM];
} Body;

void readdata(int *N, int *K, double *DT, Body bds[]) {
	int i;

	scanf("%d", N);
	scanf("%d", K);
	scanf("%lf", DT);
	for (i = 0; i < *N; ++ i) scanf("%lf", &bds[i].mass);
	for (i = 0; i < *N; ++ i) scanf("%lf %lf %lf", &(bds[i].x)[0], &(bds[i].x)[1], &(bds[i].x)[2]);
	for (i = 0; i < *N; ++ i) scanf("%lf %lf %lf", &(bds[i].v)[0], &(bds[i].v)[1], &(bds[i].v)[2]);
}

Body get_centroid(int n_body, Body bodies[]) {
	int i, k;
	Body centroid;
	centroid.mass = .0;
	for (k = 0; k < N_DIM; ++ k) centroid.x[k] = centroid.v[k] = .0;
	for (i = 0; i < n_body; ++ i) {
		centroid.mass += bodies[i].mass;
		for (k = 0; k < N_DIM; ++ k) {
			centroid.x[k] += bodies[i].mass * bodies[i].x[k];
			centroid.v[k] += bodies[i].v[k];
		}
	}
	for (k = 0; k < N_DIM; ++ k) {
		centroid.x[k] /= centroid.mass;
		centroid.v[k] /= n_body;
		centroid.mass /= n_body;
	}
	return centroid;
}

void compute_force_as_vector(const Body *this, const Body *that, double vec[]) {  // no init happens inside
	double r2, tmp;
	int k;
	
	if (this == that) return;
	r2 = pow(this->x[0] - that->x[0], 2) + pow(this->x[1] - that->x[1], 2) + pow(this->x[2] - that->x[2], 2);
	if (r2 > CUTOFF_SQR) return;
	tmp = (G * this->mass * that->mass) / pow(r2, 1.5);  // quantative variable
	for (k = 0; k < N_DIM; ++ k) vec[k] += tmp * (that->x[k] - this->x[k]);  // vector variable
}

void report(int t, double DT, const Body* c) {
	printf("\nConditions after timestep %d (time = %lf) :\n\n", t, t * DT);
	printf("\tCenter of Mass:\t(%lf, %lf, %lf)\n", c->x[0], c->x[1], c->x[2]);
	printf("\tAverage Velocity:\t(%lf, %lf, %lf)\n", c->v[0], c->v[1], c->v[2]);
}

int main(int argc, char **argv) {
	int N;  // number of bodies
	int K;  // number of time steps
	double DT;  // time step size
	
	double F[N_DIM];  // force in each dimension
	int i, j, k, t;

	Body bodies[MAX_N_BODY];
	Body centroid;

	readdata(&N, &K, &DT, bodies);
	// printf("%lf %lf %lf\n", bodies[0].x[0], bodies[0].x[1], bodies[0].x[2]);

	for (t = 0; t < K; ++ t) {
		if (t % TIME_INTERVAL == 0) {
			centroid = get_centroid(N, bodies);
			report(t, DT, &centroid);
		}
		for (i = 0; i < N; ++ i) {  // for each body
			F[0] = F[1] = F[2] = .0;
			for (j = 0; j < N; ++ j) 
				compute_force_as_vector(&bodies[i], &bodies[j], F);  // accumulates vec variable F
			for (k = 0; k < N_DIM; ++ k) {
				bodies[i].v[k] += F[k] * (DT / 2.0) / bodies[i].mass;
				bodies[i].x[k] += bodies[i].v[k] * DT; 
			}
		}
	}
	centroid = get_centroid(N, bodies);
	report(t, DT, &centroid);

	return 0;
}
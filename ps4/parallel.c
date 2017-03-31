#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include <mpi.h>

#define MAX_N_BODY 50000
#define TIME_INTERVAL 128
#define G 1.0
#define N_DIM 3
#define NP 8

typedef struct body_t {
	double mass;
	double x[N_DIM];
	double v[N_DIM];
} Body;

void sign_of_dims(int rank, int dims[]) {
	dims[2] = (rank / 4) ? 1 : -1;
	dims[1] = ((rank % 4) / 2) ? 1: -1;
	dims[0] = (rank % 2) ? : 1 : -1;
}

int get_octant_rank(int dims[]) {
	return (dims[2] > 0) * 4 + (dims[1] > 0) * 2 + (dims[0] > 0);
}

void readdata(char filename[], int *N, int *K, double *DT, Body bodies[], int send_count[], int send_displ[]) {
//		double *mass, double (*x)[N_DIM], double (*v)[N_DIM]) {
	FILE *fp;
	int i, sum, i_oct;
	Body bds[MAX_N_BODY];
	int ptrs[NP];

	if ((fp = fopen(filename, "r")) == NULL) {
		printf("Data cannot be read!\n");
		exit(1);
	}
	fscanf(fp, "%d", N);
	fscanf(fp, "%d", K);
	fscanf(fp, "%lf", DT);
	for (i = 0; i < *N; ++ i) fscanf(fp, "%lf", &bds[i].mass);
	memset(send_count, 0, NP * sizeof(int));
	for (i = 0; i < *N; ++ i) {
		fscanf(fp, "%lf %lf %lf", &(bds[i].x)[0], &(bds[i].x)[1], &(bds[i].x)[2]);
		i_oct = get_octant_rank(bds[i].x);
		send_count[i_oct] += 1;
	}
	for (i = 0; i < *N; ++ i) fscanf(fp, "%lf %lf %lf", &(bds[i].v)[0], &(bds[i].v)[1], &(bds[i].v)[2]);
	fclose(fp);
	
	sum = 0;
	for (i = 0; i < NP; ++ i) {
		send_displ[i] = sum;
		sum += send_count[i];
	}  // offset in form of 'num of Body'
	memcpy(ptrs, send_displ, NP * sizeof(int));
	
	// put each body into its octant's slot
	for (i = 0; i < *N; ++ i) {
		i_oct = get_octant_rank(bds[i].x);
		memcpy(&bodies[ptrs[i_oct]], &bds[i], sizeof(Body));
		ptrs[i_oct] += 1;
	}
}

void report(int t, int N, double DT, double mass[], double (*x)[N_DIM], double (*v)[N_DIM]) {
	double cx[N_DIM], cv[N_DIM];
	double msum = .0;
	int i, k;

	for (k = 0; k < N_DIM; ++ k) cx[k] = cv[k] = .0;
	for (i = 0; i < N; ++ i) {
		msum += mass[i];
		for (k = 0; k < N_DIM; ++ k) {
			cx[k] += mass[i] * x[i][k];
			cv[k] += v[i][k];
		}
	}
	for (k = 0; k < N_DIM; ++ k) {
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
	

	double a[N_DIM];  // acceleration in each dimension
	double F[N_DIM];  // force in each dimension
	int i, j, k, t;
	double tmp, r2;

	int octant[MAX_N_BODY];
	int dims[NP][3];
	// int n_bodies[NP];

	MPI_Status status;
	MPI_Datatype mpi_body_t;
	int rank;

	Body bodies[MAX_N_BODY];
	int n_bodies[NP];  // behave as send_count
	int send_count[NP];
	int send_displ[NP];
	Body recvbuf[MAX_N_BODY];

	MPI_Init(&argc, &argv); // Required MPI initialization call
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

	// create mpi version of Body for transmission
	MPI_Type_create_struct(3, {1, N_DIM, N_DIM}, 
		{offsetof(Body, mass), offsetof(Body, x), offsetof(Body, v)},
		{MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}, &mpi_body_t);
	MPI_Type_commit(&mpi_body_t);

	if (rank == 0) {
		readdata("data/testdata1", &N, &K, &DT, bodies, n_bodies, send_displ);
		// for (i = 0; i < NP; ++ i) {
		// 	send_count[i] = n_bodies[i] * sizeof(Body);
		// 	send_displ[i] *= sizeof(Body);
		// }
		// printf("%d %d %lf\n", N, K, DT);
		// for (i = 0; i < N; ++ i) printf("%lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);

	}

	// for (i = 0; i < 8; ++ i) sign_of_dims(i, &dims[i]);

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&DT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatterv(bodies, n_bodies, send_displ, mpi_body_t, 
		recvbuf, n_bodies[rank], mpi_body_t, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	printf("From process %d, get DT = %lf\n", rank, DT);

	// for (t = 0; t < K; ++ t) {
	// 	if (t % TIME_INTERVAL == 0) report(t, N, DT, mass, x, v);
	// 	for (i = 0; i < N; ++ i) {  // for each body
	// 		F[0] = F[1] = F[2] = .0;
	// 		for (j = 0; j < N; ++ j) {
	// 			if (j == i) continue;
	// 			r2 = pow(x[i][0] - x[j][0], 2) + pow(x[i][1] - x[j][1], 2) + pow(x[i][2] - x[j][2], 2);
	// 			if (r2 > 25.0) continue;
	// 			tmp = (G * mass[i] * mass[j]) / (r2 * sqrt(r2));  // quantative variable
	// 			for (k = 0; k < 3; ++ k) F[k] += tmp * (x[j][k] - x[i][k]);  // vector variable
	// 		}
	// 		for (k = 0; k < 3; ++ k) a[i][k] = F[k] / mass[k];
	// 		for (k = 0; k < 3; ++ k) v[i][k] += a[i][k] * (DT / 2.0);
	// 		for (k = 0; k < 3; ++ k) x[i][k] += v[i][k] * DT;
	// 	}
	// }
	// report(K, N, DT, mass, x, v);

	MPI_Type_free(&mpi_body_t);
	MPI_Finalize();
}
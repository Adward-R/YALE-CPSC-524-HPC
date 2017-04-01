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
#define CUTOFF_SQR 25.0
#define EPSILON 0.0001

typedef struct body_t {
	double mass;
	double x[N_DIM];
	double v[N_DIM];
} Body;

void sign_of_dims(int rank, int dims[]) {
	dims[2] = (rank / 4) ? 1 : -1;
	dims[1] = ((rank % 4) / 2) ? 1: -1;
	dims[0] = (rank % 2) ? 1 : -1;
}

int get_octant_rank(double dims[]) {	
	int ans = (dims[2] > 0) * 4 + (dims[1] > 0) * 2 + (dims[0] > 0);
	// printf("DIM : %lf, %lf, %lf; ANS: %d\n", dims[0], dims[1], dims[2], ans);
	return ans;
}

void readdata(int *N, int *K, double *DT, Body bodies[], int send_count[], int send_displ[]) {
	int i, sum, i_oct;
	Body bds[MAX_N_BODY];
	int ptrs[NP];

	scanf("%d", N);
	scanf("%d", K);
	scanf("%lf", DT);
	for (i = 0; i < *N; ++ i) scanf("%lf", &bds[i].mass);
	memset(send_count, 0, NP * sizeof(int));
	for (i = 0; i < *N; ++ i) {
		scanf("%lf %lf %lf", &(bds[i].x)[0], &(bds[i].x)[1], &(bds[i].x)[2]);
		i_oct = get_octant_rank(bds[i].x);
		send_count[i_oct] += 1;
	}
	for (i = 0; i < *N; ++ i) scanf("%lf %lf %lf", &(bds[i].v)[0], &(bds[i].v)[1], &(bds[i].v)[2]);
	
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

void report(int t, double DT, const Body* c, int n_bodies[]) {
	int i;
	// printf("\n\nFrom process %d...", rank);
	printf("\nConditions after timestep %d (time = %lf) :\n\n", t, t * DT);
	printf("\tCenter of Mass:\t(%lf, %lf, %lf)\n", c->x[0], c->x[1], c->x[2]);
	printf("\tAverage Velocity:\t(%lf, %lf, %lf)\n", c->v[0], c->v[1], c->v[2]);
	printf("\tBodies in octants' distribution:");
	for (i = 0; i < NP; ++ i) printf(" %d", n_bodies[i]);
}

int octant_dims[][N_DIM] = {
	{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
	{-1, -1, 1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 1}
};

void shouldSend(double pos[N_DIM], int transmit[NP]) {
	int o, k;
	double dist;

	for (o = 0; o < NP; ++ o) {
		dist = .0;
		for (k = 0; k < N_DIM; ++ k)
			if (octant_dims[o][k] * pos[k] < .0) dist += pow(pos[k], 2);
		if (abs(dist) < EPSILON || dist > CUTOFF_SQR) transmit[o] = 0;
		else transmit[o] = 1;
	}
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


int main(int argc, char **argv) {
	int N;  // number of bodies
	int K;  // number of time steps
	double DT;  // time step size
	
	double a[N_DIM];  // acceleration in each dimension
	double F[N_DIM];  // force in each dimension
	int i, j, k, t, cnt;
	double sum;

	MPI_Status status;
	MPI_Datatype mpi_body_t;
	int rank, i_oct;

	Body bodies[MAX_N_BODY];
	int n_bodies[NP];  // behave as send_count
	int send_count[NP];
	int send_displ[NP];
	int transmit[NP];
	int recv_count[NP];
	int recv_displ[NP];
	Body sendbuf[MAX_N_BODY];
	Body recvbuf[MAX_N_BODY];
	int ptrs[NP];
	Body centroid;

	int block_length_array[] = {1, N_DIM, N_DIM};
	MPI_Aint displ_array[] = {offsetof(Body, mass), offsetof(Body, x), offsetof(Body, v)};
	MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Init(&argc, &argv); // Required MPI initialization call
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

	// create mpi version of Body for transmission
	MPI_Type_create_struct(3, block_length_array, displ_array, types, &mpi_body_t);
	MPI_Type_commit(&mpi_body_t);

	
	readdata(&N, &K, &DT, bodies, n_bodies, send_displ);

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&DT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(n_bodies, NP, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// for (i = 0; i < NP; ++ i) printf("%d ", n_bodies[i]);
	// printf("\n");
	MPI_Scatterv(bodies, n_bodies, send_displ, mpi_body_t, 
		recvbuf, n_bodies[rank], mpi_body_t, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	memcpy(bodies, recvbuf, n_bodies[rank] * sizeof(Body));  // release recvbuf
	

	// printf("\n---Proc %d---\n", rank);
	// 	for (i = 0; i < n_bodies[rank]; ++ i) {
	// 		for (j = 0; j < 3; ++ j) printf("%lf ", bodies[i].x[j]);
	// 		printf("\n----\n");
	// 	}
	

	for (t = 0; t < K; ++ t) {
		if (t % TIME_INTERVAL == 0) {
		// if (t % TIME_INTERVAL == 0) {
			// printf("%lf %lf %lf\n", bodies[0].x[0], bodies[0].x[1], bodies[0].x[2]);
			centroid = get_centroid(n_bodies[rank], bodies);
			// printf("%lf %lf %lf\n", centroid.x[0], centroid.x[1], centroid.x[2]);
			MPI_Gather(&centroid, 1, mpi_body_t, recvbuf, 1, mpi_body_t, 0, MPI_COMM_WORLD);
			if (rank == 0) {
				centroid = get_centroid(NP, recvbuf);
				report(t, DT, &centroid, n_bodies);
			}
		}

		/* Prepare sending bodies within 5DU to other octants */
		memset(send_count, 0, NP * sizeof(int));

		// Compute send_count & set slots for send buffer fulfilling
		for (i = 0; i < n_bodies[rank]; ++ i) {
			shouldSend(bodies[i].x, transmit);
			for (j = 0; j < NP; ++ j) {
				if (transmit[j]) send_count[j] += 1;
			}
		}

		sum = 0;
		for (j = 0; j < NP; ++ j) {
			send_displ[j] = sum;
			sum += send_count[j];
		}
		memcpy(ptrs, send_displ, NP * sizeof(int));
		// Send buffer fulfilling
		for (i = 0; i < n_bodies[rank]; ++ i) {
			shouldSend(bodies[i].x, transmit);
			for (j = 0; j < NP; ++ j) {
				if (transmit[j]) {
					memcpy(&sendbuf[ptrs[j]], &bodies[i], sizeof(Body));
					ptrs[j] += 1;
				}
			}
		}

		// Notice each other about send_count as recv_count
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// Constitutes recv_displ
		sum = 0;
		for (j = 0; j < NP; ++ j) {
			recv_displ[j] = sum;
			sum += recv_count[j];
		}

		// Sending actual bodies
		MPI_Alltoallv(sendbuf, send_count, send_displ, mpi_body_t,
			recvbuf, recv_count, recv_displ, mpi_body_t, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		memset(send_count, 0, NP * sizeof(int));
		for (i = 0; i < n_bodies[rank]; ++ i) {  // for each body
			F[0] = F[1] = F[2] = .0;
			for (j = 0; j < n_bodies[rank]; ++ j) 
				compute_force_as_vector(&bodies[i], &bodies[j], F);  // accumulates vec variable F
			for (j = 0; j < recv_displ[NP-1] + recv_count[NP-1]; ++ j)
				compute_force_as_vector(&bodies[i], &recvbuf[j], F);
			for (k = 0; k < N_DIM; ++ k) {
				bodies[i].v[k] += F[k] * (DT / 2.0) / bodies[i].mass;
				bodies[i].x[k] += bodies[i].v[k] * DT; 
			}

			// check if the body is still within reign of current octant
			i_oct = get_octant_rank(bodies[i].x);
			if (i_oct != rank) send_count[i_oct] += 1;  // get size of slots for sendbuf
		}

		/* Communicate with each other to exchange bodies beyond their own reigns */
		// Prepare slots for sendbuf
		sum = 0;
		for (j = 0; j < NP; ++ j) {
			send_displ[j] = sum;
			sum += send_count[j];
		}
		memcpy(ptrs, send_displ, NP * sizeof(int));
		// Send buffer fulfilling
		cnt = n_bodies[rank];  // record original belonging number
		for (i = 0; i < cnt; ++ i) {
			i_oct = get_octant_rank(bodies[i].x);
			if (rank == i_oct) continue;
			memcpy(&sendbuf[ptrs[i_oct]], &bodies[i], sizeof(Body));
			ptrs[i_oct] += 1;
			// bodies[i].mass = FP_NAN;
			// n_bodies[rank] -= 1;
		}

		// Notice each other about send_count as recv_count
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// Constitutes recv_displ
		sum = 0;
		for (j = 0; j < NP; ++ j) {
			recv_displ[j] = sum;
			sum += recv_count[j];
		}

		// printf("\n---\nIn process %d~", rank);
		// printf("\nSend count: ");
		// for (i = 0; i < NP; ++ i) printf("%d ", send_count[i]);
		// printf("\nRecv count: ");
		// for (i = 0; i < NP; ++ i) printf("%d ", recv_count[i]);
		// printf("---\n");

		// Sending actual bodies
		MPI_Alltoallv(sendbuf, send_count, send_displ, mpi_body_t,
			recvbuf, recv_count, recv_displ, mpi_body_t, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// put each received body into missing slots
		i = recv_displ[NP-1] + recv_count[NP-1];
		for (j = 0; j < cnt; ++ j) {
			if (rank == get_octant_rank(bodies[j].x)) {
				memcpy(&recvbuf[i], &bodies[j], sizeof(Body));
				i += 1;
			}
		}
		memcpy(bodies, recvbuf, i * sizeof(Body));
		n_bodies[rank] = i;

		// printf("%lf\n", bodies[0].x[0]);
		
		// i = j = 0;
		// for (; i < recv_displ[NP-1] + recv_count[NP-1]; ++ i) {
		// 	while (j < cnt && !isnan(bodies[j].mass)) j += 1;
		// 	memcpy(&bodies[j], &recvbuf[i], sizeof(Body));
		// 	n_bodies[rank] += 1;
		// 	j += 1;
		// }
		// // n_bodies[rank] += recv_displ[NP-1] + recv_count[NP-1];
		// if (j < cnt) {
		// 	for (i = j; i < cnt; ++ i) {
		// 		if (!isnan(bodies[i].mass)) {
		// 			memcpy(&bodies[j], &bodies[i], sizeof(Body));
		// 			j += 1;
		// 		}
		// 	}
		// }

	}

	centroid = get_centroid(n_bodies[rank], bodies);
	MPI_Gather(&centroid, 1, mpi_body_t, recvbuf, 1, mpi_body_t, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		centroid = get_centroid(NP, recvbuf);
		report(t, DT, &centroid, n_bodies);
	}

	MPI_Type_free(&mpi_body_t);
	MPI_Finalize();
}
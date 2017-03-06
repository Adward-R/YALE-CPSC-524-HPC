#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include <mpi.h>

#define CACHE_SIZE 64000000
#define BLK_SIZE 20000000
#define TAG_ROW 11
#define TAG_COL 22
#define TAG_RESULT 33
#define MIN(a,b) (((a)<(b))?(a):(b))

int prev(int col, int limit) { return (col-1>=0) ? col-1 : limit-1; }
int succ(int col, int limit) { return (col+1<limit) ? col+1 : 0; }
       
int group_length(int size, int idx) {
  /* row-group i contains row [i*size, (i+1)*size), row j in lower triangle A has j+1 cells,
    so row-group i contains: size * (i*size + (1+size)/2) cells. */
  return idx * size * size + (1+size) * size / 2;
}

int group_index_start(int size, int idx) {
  /* group [0,i) contains overall: start_index_group(i) = (i*size/2) * (i*size + 1) cells,
    group i's index range: [start_index_group(i), start_index_group(i) + group_size(i)). */
  return (idx * size / 2) * (idx * size + 1);
}

double dense_matmul(int N, int size, int gi, int gj, double gA[], double gB[], double ret[]) {
  int i, j, si, sj, k, dst, src;
  double gC[BLK_SIZE];
  double wct;
  // gA(k) countains row [k*size, (k+1)*size)
  wct = MPI_Wtime();
  for (i = 0; i < size; ++i) {
    for (j = 0; j < size; ++j) {
      gC[i*size+j] = 0.0;
      si = gi * size * i + (1 + i) * i / 2;
      sj = gj * size * j + (1 + j) * j / 2;
      for (k = 0; k < MIN(gi * size + i, gj * size + j) + 1; ++k) {
        gC[i*size+j] += gA[si+k] * gB[sj+k];
      }
    }
  }

  for (i = 0; i < size; ++i) {
      dst = i * N + gj * size; // gj <- col in the outer scope
      src = i * size;
      for (k = 0; k < size; ++k) { ret[dst+k] = gC[src+k]; }
  }
  return MPI_Wtime() - wct;
}

int main(int argc, char **argv) {

  int N, i, j, k, run, si, sj;
  double *A, *B, *C;
  int sizeAB, sizeC, dst, src;
  double cpt = 0.0, cmt = 0.0, wct; // compute-time, communicate-time, wall-clock-time
  int col; // col-group index

  MPI_Status status;

  int np, size, rank;
  int g_start, g_len;

  double gA[CACHE_SIZE], gB[CACHE_SIZE], ret[CACHE_SIZE];
  double nB[CACHE_SIZE]; // used as recving cache of column groups

  if (argc <= 1) {
    printf("Please specify size of multiplied matrices!\n");
    exit(1);
  }
  N = atoi(argv[1]);

  MPI_Init(&argc, &argv); // Required MPI initialization call
  MPI_Comm_size(MPI_COMM_WORLD, &np); // Get no. of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

  sizeAB = N * (N + 1) / 2; //Only enough space for the nonzero portions of the matrices
  sizeC = N * N; // All of C will be nonzero, in general!
  size = N / np; // width of grouped rows or cols

  if (rank == 0) {
    if (N == 1000) {
      printf("Matrix multiplication times (secs):\nN\tCompute\tCommunicate\n");
      printf("-----\t-------------\t-------------\n");
    }
    A = (double *) calloc(sizeAB, sizeof(double));
    B = (double *) calloc(sizeAB, sizeof(double));
    C = (double *) calloc(sizeC, sizeof(double));

    srand(12345); // Use a standard seed value for reproducibility
    // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
    for (i=0; i<sizeAB; i++) A[i] = ((double) rand() / (double) RAND_MAX);
    for (i=0; i<sizeAB; i++) B[i] = ((double) rand() / (double) RAND_MAX);
  
    for (i = 0; i < np; ++i) {
      g_start = group_index_start(size, i);
      g_len = group_length(size, i);
      for (j = 0; j < g_len; ++j) { gA[j] = A[g_start + j]; }
      for (j = 0; j < g_len; ++j) { gB[j] = B[g_start + j]; }
      wct = MPI_Wtime();
      MPI_Send(gA, g_len, MPI_DOUBLE, i, TAG_ROW, MPI_COMM_WORLD);
      MPI_Send(gB, g_len, MPI_DOUBLE, i, TAG_COL, MPI_COMM_WORLD);
      cmt += MPI_Wtime() - wct;
    }
    
  }
  
  g_len = group_length(size, rank);
  wct = MPI_Wtime();
  MPI_Recv(gA, g_len, MPI_DOUBLE, 0, TAG_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(gB, group_length(size, rank), MPI_DOUBLE, 0, TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  cmt += MPI_Wtime() - wct;

  col = rank;
  for (j = 0; j < np; ++j) {
    cpt += dense_matmul(N, size, rank, col, gA, gB, ret);
    wct = MPI_Wtime();
    if (rank % 2) {
      MPI_Send(gB, group_length(size, col), MPI_DOUBLE, succ(rank, np), TAG_COL, MPI_COMM_WORLD);
      MPI_Recv(nB, group_length(size, prev(col, np)), MPI_DOUBLE, prev(rank, np), TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(nB, group_length(size, prev(col, np)), MPI_DOUBLE, prev(rank, np), TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(gB, group_length(size, col), MPI_DOUBLE, succ(rank, np), TAG_COL, MPI_COMM_WORLD);
    }
    cmt += MPI_Wtime() - wct;

    col = prev(col, np);
    memcpy(gB, nB, CACHE_SIZE * sizeof(double));
  }

  wct = MPI_Wtime();
  MPI_Send(ret, size * N, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
  cmt += MPI_Wtime() - wct;

  if (rank == 0) {
    wct = MPI_Wtime();
    for (i = 0; i < np; ++i) {
      MPI_Recv(ret, size * N, MPI_DOUBLE, i, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      dst = i * size * N;
      src = 0;
      for (k = 0; k < size * N; ++k) { C[dst+k] = ret[src+k]; }
    }
    cmt += MPI_Wtime() - wct;

    printf("%5d\t%9.4lf\t%9.4lf\n", N, cpt, cmt);

    free(A);
    free(B);
    free(C);
  }

  MPI_Finalize();
}

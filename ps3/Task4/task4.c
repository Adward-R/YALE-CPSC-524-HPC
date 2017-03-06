#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include <mpi.h>

#define CACHE_SIZE 64000000 // max((N*N)/(2*NP)) 32000000
#define BLK_SIZE 20000000 // 2000000
#define TAG_ROW 11
#define TAG_COL 22
#define TAG_RESULT 33
#define TAG_PARTITION 44
#define MAX_NP 8
#define MAX_SIZE 4000
#define MIN(a,b) (((a)<(b))?(a):(b))

int prev(int col, int limit) { return (col-1>=0) ? col-1 : limit-1; }
int succ(int col, int limit) { return (col+1<limit) ? col+1 : 0; }
int rc_start(int rc_num) { return rc_num * (rc_num + 1) / 2; }

double dense_matmul(int N, int psize_row, int psize_col, 
    int row_partition[], int col_partition[], double gA[], double gB[], double ret[]) {
  int i, j, si, sj, k, dst, src;
  double gC[BLK_SIZE];
  double wct;

  wct = MPI_Wtime();
  dst = 0;
  si = 0;
  for (i = 0; i < psize_row; ++i) {
    sj = 0;
    for (j = 0; j < psize_col; ++j) {
      gC[dst] = 0.0;
      for (k = 0; k <= MIN(row_partition[i], col_partition[j]); ++k) {
        gC[dst] += gA[si+k] * gB[sj+k];
      }
      dst += 1;
      sj += col_partition[j] + 1;
    }
    si += row_partition[i] + 1;
  }
  
  src = 0;
  for (i = 0; i < psize_row; ++i) {
    for (j = 0; j < psize_col; ++j) {
      ret[i*N + col_partition[j]] = gC[src++];
    }
  }

  return MPI_Wtime() - wct;
}

double MPI_Wait_Timed(MPI_Request *request) {
  double wct;
  wct = MPI_Wtime();
  MPI_Wait(request, MPI_STATUS_IGNORE);
  return MPI_Wtime() - wct;
}

int main(int argc, char **argv) {

  int N, i, j, k, si, sj, ip, sum;
  double *A, *B, *C;
  int sizeAB, sizeC, dst, src, flag, counter;
  double cpt = 0.0, cmt = 0.0, wct; // compute-time, communicate-time, wall-clock-time
  int col; // col-group index

  MPI_Status status;
  MPI_Request send_reqs[MAX_NP][2], recv_reqs[MAX_NP][2]; // [0] for data, [1] for partition number

  int np, rank;
  int size, length, last_size, last_length;
  int rank_len, g_len, p_len, prev_gl, prev_pl;

  double gA[CACHE_SIZE], gB[CACHE_SIZE], ret[CACHE_SIZE];
  double nB[CACHE_SIZE]; // used as recving cache of column groups
  int partitions[MAX_NP][MAX_SIZE];
  int row_partition[MAX_SIZE], col_partition[MAX_SIZE], new_col_part[MAX_SIZE];
  int psizes[MAX_NP];

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
  size = (int) ((N - 1) / np) + 1; // width of grouped rows or cols
  length = (int) ((N * N - 1) / np) + 1;
  last_size = N - size * (np-1);
  last_length = N * N - length * (np - 1);

  // printf("last_length: %d\n", last_length);
  // printf("last_size: %d\n", last_size);

  if (rank == 0) {
    if (N == 4000) {
      printf("Matrix multiplication times (secs):\nN\tCompute\t     Communicate\n");
      printf("-----\t-------------\t-------------\n");
    }
    A = (double *) calloc(sizeAB, sizeof(double));
    B = (double *) calloc(sizeAB, sizeof(double));
    C = (double *) calloc(sizeC, sizeof(double));

    srand(12345); // Use a standard seed value for reproducibility
    // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
    for (i=0; i<sizeAB; i++) A[i] = ((double) rand() / (double) RAND_MAX);
    for (i=0; i<sizeAB; i++) B[i] = ((double) rand() / (double) RAND_MAX);
  
    i = 0;
    j = N - 1;
    for (ip = 0; ip < np; ++ip) {
      k = 0; // partition cache pointer
      while (i <= j) {
        if (i <= j && k < size) {
          partitions[ip][k] = i;
          k += 1;
          i += 1;
        } else break;
        if (i <= j && k < size) {
          partitions[ip][k] = j;
          k += 1;
          j -= 1;
        } else break;
      }
      psizes[ip] = k;
      // printf("psizes: %d\n", psizes[ip]);
    }

    // copy assigned rows and cols to sending caches
    for (ip = 0; ip < np; ++ip) {
      dst = 0;
      for (k = 0; k < psizes[ip]; ++k) { // copy row by row & col by col
        src = rc_start(partitions[ip][k]);
        for (i = 0; i <= partitions[ip][k]; ++i) {
          gA[dst + i] = A[src + i];
          gB[dst + i] = B[src + i];
        }
        dst += partitions[ip][k] + 1;
      }
      wct = MPI_Wtime();
      MPI_Send(gA, dst, MPI_DOUBLE, ip, TAG_ROW, MPI_COMM_WORLD);
      MPI_Send(gB, dst, MPI_DOUBLE, ip, TAG_COL, MPI_COMM_WORLD);
      MPI_Send(partitions[ip], psizes[ip], MPI_DOUBLE, ip, TAG_PARTITION, MPI_COMM_WORLD);
      cmt += MPI_Wtime() - wct;
    }
    
  } // end of master process, part-I: task distribution

  g_len = (rank < np-1) ? length : last_length;
  p_len = (rank < np-1) ? size : last_size;
  rank_len = p_len;
  wct = MPI_Wtime();
  MPI_Recv(gA, g_len, MPI_DOUBLE, 0, TAG_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(gB, g_len, MPI_DOUBLE, 0, TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(row_partition, p_len, MPI_DOUBLE, 0, TAG_PARTITION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  memcpy(col_partition, row_partition, p_len * sizeof(double));
  cmt += MPI_Wtime() - wct;

  col = rank;
  for (j = 0; j < np; ++j) {
    prev_gl = (col != 0) ? length : last_length;
    prev_pl = (col != 0) ? size : last_size;
    MPI_Irecv(nB, prev_gl, MPI_DOUBLE, prev(rank, np), TAG_COL, MPI_COMM_WORLD, &recv_reqs[rank][0]);
    MPI_Irecv(new_col_part, prev_pl, MPI_DOUBLE, prev(rank, np), TAG_PARTITION, MPI_COMM_WORLD, &recv_reqs[rank][1]);

    // doing work using old data, while receiving goes parallel
    p_len = (col < np-1) ? size : last_size;
    cpt += dense_matmul(N, rank_len, p_len, row_partition, col_partition, gA, gB, ret);

    g_len = (col < np-1) ? length : last_length;
    // wait for sending to complete so that queue cache will not burst
    // printf("wkr %d sending to %d\n", rank, succ(rank, np));
    MPI_Isend(gB, g_len, MPI_DOUBLE, succ(rank, np), TAG_COL, MPI_COMM_WORLD, &send_reqs[rank][0]);
    cmt += MPI_Wait_Timed(&send_reqs[rank][0]);
    MPI_Isend(col_partition, p_len, MPI_DOUBLE, succ(rank, np), TAG_PARTITION, MPI_COMM_WORLD, &send_reqs[rank][1]);    
    cmt += MPI_Wait_Timed(&send_reqs[rank][1]);

    col = prev(col, np);
    cmt += MPI_Wait_Timed(&recv_reqs[rank][0]);
    memcpy(gB, nB, prev_gl * sizeof(double));
    
    cmt += MPI_Wait_Timed(&recv_reqs[rank][1]);
    memcpy(col_partition, new_col_part, prev_pl * sizeof(double));
  }
  
  // merely sending back data results is sufficient, cuz master knows orginal row partition
  p_len = (rank < np-1) ? size : last_size;
  MPI_Send(ret, p_len * N, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);

  if (rank == 0) {
    wct = MPI_Wtime();
    for (ip = 0; ip < np; ++ip) {
      p_len = (ip < np-1) ? size : last_size;
      MPI_Recv(ret, p_len * N, MPI_DOUBLE, ip, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (i = 0; i < p_len; ++i) {
        dst = partitions[ip][i] * N;
        src = i * N;
        for (k = 0; k < N; ++k) { C[dst+k] = ret[src+k]; }
      }
    }
    cmt += MPI_Wtime() - wct;

    printf("%5d\t%9.4lf\t%9.4lf\n", N, cpt, cmt);

    // for (k = 0; k < N; k++) { printf("%lf\n", C[k*N+k]);}

    free(A);
    free(B);
    free(C);
  }

  MPI_Finalize();
}

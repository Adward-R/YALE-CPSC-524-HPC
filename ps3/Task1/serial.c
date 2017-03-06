#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
double matmul(int, double*, double*, double*);

int main(int argc, char **argv) {

  /*
    This is the serial main program for CPSC424/524 Assignment #3.

    Author: Andrew Sherman, Yale University

    Date: 2/01/2016

  */

  int N, i, j, k, run;
  double *A, *B, *C;
  int sizeAB, sizeC, iA, iB, iC;

  int sizes[5]={1000,2000,4000,8000,12000};

  double wctime;

  printf("Matrix multiplication times:\n   N      TIME (secs)\n -----   -------------\n");
  for (run=0; run<5; run++) {
    N = sizes[run];

    sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrices
    sizeC = N*N; // All of C will be nonzero, in general!

    A = (double *) calloc(sizeAB, sizeof(double));
    B = (double *) calloc(sizeAB, sizeof(double));
    C = (double *) calloc(sizeC, sizeof(double));
  
    srand(12345); // Use a standard seed value for reproducibility

    // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
    for (i=0; i<sizeAB; i++) A[i] = ((double) rand()/(double)RAND_MAX);
    for (i=0; i<sizeAB; i++) B[i] = ((double) rand()/(double)RAND_MAX);

    wctime = matmul(N, A, B, C);

    printf ("  %5d    %9.4f\n", N, wctime);

    // for (i = 0; i < N; ++i) {
    //   for (j = 0; j < N; ++j) {
    //     printf("%d ", C[i*N+j]);
    //   }
    //   printf("\n");
    // }

    free(A);
    free(B);
    free(C);
  }

}

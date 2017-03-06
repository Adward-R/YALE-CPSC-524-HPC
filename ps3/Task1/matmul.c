#include "timing.h"
#define MIN(a,b) (((a)<(b))?(a):(b))

double matmul(int N, double* A, double* B, double* C) {

/*
  This is the serial version of triangular matrix multiplication for CPSC424/524 Assignment #3.

  Author: Andrew Sherman, Yale University

  Date: 2/01/2016

*/

  int i, j, k;
  int iA, iB, iC;
  double wctime0, wctime1, cputime;

  timing(&wctime0, &cputime);

// This loop computes the matrix-matrix product
  iC = 0;
  for (i=0; i<N; i++) {
    iA = i*(i+1)/2; // initializes row pointer in A
    for (j=0; j<N; j++,iC++) {
      iB = j*(j+1)/2; // initializes column pointer in B
      C[iC] = 0.;
      for (k=0; k<=MIN(i,j); k++) C[iC] += A[iA+k] * B[iB+k]; // avoids using known-0 entries 
    }
  }

  timing(&wctime1, &cputime);
  return(wctime1 - wctime0);
}

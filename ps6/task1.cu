#define FP float

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <math.h>

__global__ void gpu_matrixmult(FP *a,FP *b, FP *c, int n, int p, int m) {

    int row = threadIdx.y + blockDim.y * blockIdx.y;
    int col = threadIdx.x + blockDim.x * blockIdx.x;

    int indexb = col;
    int index = row * m + col;

    if(col < m && row < n) {
        c[index] = 0.;
        for (int indexa = row * p; indexa < (row * p + p); indexa ++, indexb += m)
            c[index] += a[indexa] * b[indexb];
    }

}

void cpu_matrixmult(FP *a,FP *b, FP *c, int n, int p, int m) {

    int i, j, k;
    for (k = 0; k < p; k++) {
        for (i = 0; i < n; i++) {
            FP r = a[i*p+k];
            for (j = 0; j < m; j++) {
                c[i*m+j] += r * b[k*m+j];
            }
        }
    }
}

int main(int argc, char *argv[]) {

    int i, j; // loop counters

    int gpucount = 0; // Count of available GPUs
    int gpunum = 0; // Device number to use

    int n, p, m; // matrix dimension
    FP *a, *b, *c; // 1d-array representing matrices, where a(nxp), b(pxm), c(nxm).
    FP *dev_a, *dev_b, *dev_c;

    cudaEvent_t start, stop; // using cuda events to measure time
    float elapsed_time_ms; // which is applicable for asynchronous code also
    cudaError_t errorcode;

    // --------------------SET PARAMETERS AND DATA -----------------------

    errorcode = cudaGetDeviceCount(&gpucount);
    if (errorcode == cudaErrorNoDevice) {
        printf("No GPUs are visible\n");
        exit(-1);
    } else {
        printf("Device count = %d\n", gpucount);
    }

    if ((argc < 4) || (argc > 7)) {
        printf("Usage: matmul <matrix dim> <block dim> <grid dim> [<dev num>]\n");
        exit(-1);
    }

    n = atoi(argv[1]);
    p = atoi(argv[2]);
    m = atoi(argv[3]);

//    Block_Dim = 32; // Block_dim ^ 2 <= 1024
//    Grid_Dim = 1; // Grid_dim * Block_dim >= MAX(m, n)
//    int max_size = (n > m) ? n : m;
//    while (Block_Dim * Grid_Dim < max_size) Grid_Dim *= 2;
//    Grid_Dim = 1024;

    if (argc == 5) {
        gpunum = atoi(argv[4]); // Device number
        if ((gpunum > 2) || (gpunum < 0)) {
            printf("Error, Device number must be 0, 1, or 2\n");
            exit(-1);
        }
    }
    cudaSetDevice(gpunum);
    printf("Using device %d\n", gpunum);

    int Grid_dim_x, Grid_dim_y;
    int Block_dim_x, Block_dim_y;

    Block_dim_x = Block_dim_y = 32;
    Grid_dim_x = m / Block_dim_x;
    Grid_dim_y = n / Block_dim_y;

    printf("Matrix Dimension = %d x %d\n", n, m);
    printf("Block_Dim = %d x %d, Grid_Dim = %d x %d\n", Block_dim_x, Block_dim_y, Grid_dim_x, Grid_dim_y);

    dim3 Grid(Grid_dim_x, Grid_dim_y); //Grid structure
    dim3 Block(Block_dim_x, Block_dim_y); //Block structure

    a = (FP *) malloc(n * p * sizeof(FP)); // dynamically allocated memory for arrays on host
    b = (FP *) malloc(p * m * sizeof(FP));
    c = (FP *) malloc(n * m * sizeof(FP)); // results from GPU

    srand(12345);
    for (i = 0; i < n; i++)
        for (j = 0; j < p; j++) {
            a[i * p + j] = (FP) rand() / (FP) RAND_MAX;
            // a[i * p + j] = (FP) i + j; // may be helpful for debugging
        }

    for (i = 0; i < p; i++)
        for (j = 0; j < m; j++) {
            b[i * m + j] = (FP) rand() / (FP) RAND_MAX;
            // b[i * m + j] = (FP) i + j; // may be helpful for debugging
        }

    // ------------- COMPUTATION DONE ON GPU ----------------------------

    cudaMalloc((void **) &dev_a, n * p * sizeof(FP)); // allocate memory on device
    cudaMalloc((void **) &dev_b, p * m * sizeof(FP));
    cudaMalloc((void **) &dev_c, n * m * sizeof(FP));

    cudaMemcpy(dev_a, a, n * p * sizeof(FP), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, p * m * sizeof(FP), cudaMemcpyHostToDevice);

    cudaEventCreate(&start); // instrument code to measure start time
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    // cudaEventSynchronize(start); // not needed

    gpu_matrixmult << < Grid, Block >> > (dev_a, dev_b, dev_c, n, p, m);

    cudaEventRecord(stop, 0); // instrument code to measure end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    cudaMemcpy(c, dev_c, n * m * sizeof(FP), cudaMemcpyDeviceToHost);
    for (i = 0; i < n; i+=32) {
        // for (j = 0; j < m; j++) {
        printf("%f ", c[i * m + i]);
        // }
    }
    printf("\n");

    printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time

    // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
    // DEBUGGING USE ONLY (AND FOR LIMITED NUMBERS OF TIMING RUNS)

    cudaEventRecord(start, 0); // use same timing
    // cudaEventSynchronize(start); // not needed

    for (i = 0; i < n * m; ++ i) c[i] = 0;
    cpu_matrixmult(a, b, c, n, p, m); // do calculation on host (NOTE: This computes the diff with GPU result.)

    cudaEventRecord(stop, 0); // instrument code to measue end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop );

    printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms); // exec. time

// ------------------- check device creates correct results -----------------

  double error, suma, sumb, sumc, ai, bi, ci;
  suma = 0.; sumb = 0; sumc = 0;
  for(i=0;i < n*p;i++) {
    ai = (double) a[i];
    suma += ai*ai;
  }
  for(i=0;i < p*m;i++) {
    bi = (double) b[i];
    sumb += bi*bi;
  }
  for(i=0;i < n*m;i++) {
    ci = (double) c[i];
    sumc += ci*ci;
  }
  suma = sqrt(suma);
  sumb = sqrt(sumb);
  sumc = sqrt(sumc);
  error =  sumc/(n*suma*sumb);
  printf("Scaled error between GPU and CPU: %e\n", error);

// -------------- clean up ---------------------------------------

    free(a);
    free(b);
    free(c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}

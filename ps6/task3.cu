#define FP float
#define TW 32
#define TH 7

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <math.h>

__global__ void gpu_matmul_tiled(FP *a,FP *b, FP *c, int n, int p, int m) {

    FP cvalues[TH];
    __shared__ FP atile[TW][TW], btile[TH][TW][TW];

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int cols[TH]; // int col = tx + blockDim.x * blockIdx.x;
    int row = ty + blockDim.y * blockIdx.y;

    // init
    for (int i = 0; i < TH; ++i ) {
        cols[i] = tx + blockDim.x * (TH * blockIdx.x + i);
        cvalues[i] = 0;
    }

    // loop over tiles
    for (int i = 0; i < p / TW; ++ i) {
        atile[ty][tx] = a[row*p + i*TW + tx]; //Copy to shared memory
        for (int j = 0; j < TH; ++ j) {
            if (cols[j] < m) {  //Copy to shared memory
                btile[j][ty][tx] = b[(i * TW + ty) * m + cols[j]];
            }
        }
        __syncthreads();

        for (int j = 0; j < TH; ++ j) {
            for (int k = 0; k < TW; ++ k) {
                cvalues[j] += atile[ty][k] * btile[j][k][tx];
            }
        }
        __syncthreads();
    }

    // copy back to shared results
    for (int j = 0; j < TH; ++ j) {
        if (cols[j] < m) {
            c[row * m + cols[j]] = cvalues[j];
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

    Block_dim_x = Block_dim_y = TW;
    Grid_dim_x = m / Block_dim_x / TH;
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

    // cudaFuncSetCacheConfig(gpu_matmul_tiled, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(gpu_matmul_tiled, cudaFuncCachePreferShared);

    gpu_matmul_tiled << < Grid, Block >> > (dev_a, dev_b, dev_c, n, p, m);

    cudaEventRecord(stop, 0); // instrument code to measure end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    cudaMemcpy(c, dev_c, n * m * sizeof(FP), cudaMemcpyDeviceToHost);

    // Printing out diagonal to validate correctness
    for (i = 0; i < n; i+=32) {
        printf("%f ", c[i * m + i]);
        // }
    }
    printf("\n");

    printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time

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

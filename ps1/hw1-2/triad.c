#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern void timing(double* wcTime, double* cpuTime);
extern void dummy();

void kernel(int N, double a[], double b[], double c[], double d[]) {
    int i = 0;
    for (; i < N; i++) {
        a[i] = b[i] + c[i] * d[i];
    }
}

double* init(int N) {
    double* arr = (double *) malloc(sizeof(double) * N);
    double drand_max = 100.0 / (double) RAND_MAX;
    int i = 0;
    for (; i < N; i++) {
        arr[i] = drand_max * (double) rand();
    }
    return arr;
}

int main() {
    int k, r, repeat;
    long N;
    double wcs, wce, ct, runtime;
    double *a, *b, *c, *d;

    for (k = 3; k <= 24; k++) {
        N = floor(pow(2.1, k));
        a = init(N), b = init(N), c = init(N), d = init(N);
        runtime = 0.;
        repeat = 1;
        while (runtime < 1.0) {
            timing(&wcs, &ct);
            for (r = 0; r < repeat; r++) {
                kernel(N, a, b, c, d);
                if (a[N>>1] < 0.) { dummy(); }  // fool the compiler
            }
            timing(&wce, &ct);
            runtime = wce - wcs;
            repeat *= 2;
        }
        repeat /= 2;
        // 2 * N * repeat FP operations in total
        printf("N=%ld, log(N) = %d, MFlops/s = %2.1lf G, repeat=%d, runtime=%lf\n", N, k, (2 * N * repeat/pow(10, 9)) / runtime, repeat, runtime);

    }
    return 0;
}

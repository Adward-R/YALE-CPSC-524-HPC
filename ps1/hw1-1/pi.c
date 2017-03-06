#include <stdio.h>
#include <math.h>
#define N 1000000000
extern void timing(double* wcTime, double* cpuTime);

double approx_pi() {
    double sum = 0., slice = 1.0 / N, x;
    int i;
    x = slice / 2;
    for (i = 1; i < N; i++) {  // each cycle contains 6 FP operations
        sum += slice / (1.0 + x * x);
        x += slice;
    }
    return 4.0 * sum; // total FP operations: 4N
}

int main() {
    double t1, t2, cpu_t, pi;
    timing(&t1, &cpu_t);
    pi = approx_pi();
    printf("Pi=%lf, sin(Pi)=%lf\n", pi, sin(pi));
    timing(&t2, &cpu_t);
    printf("Time spent: %lf, ", t2 - t1);
    printf("MFlops/s = %lf\n", (4 * (N / 1000000)) / (t2 - t1));
    return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> //gettimeofday()
#include <time.h>
#include <omp.h>

// gcc dp_so.c -o dp_so -fopenmp
#define DATA_TYPE float

const int N = 1e9;

int main ()
{
    int i, nthreads, tid;
    DATA_TYPE x_seq, x_par, *y, *z;
    struct timeval time;
    double tstart_cpu, tend_cpu, tstart_wall, tend_wall;
    double walltime_seq, walltime_par, cputime_seq, cputime_par;

    nthreads = 8;
    printf("- - -DOT PROCUCT: OPENMP - - -\n");
    printf("Vector size           : %d\n", N);
    printf("Number of threads used: %d\n", nthreads);


    // INITIALIZATION

    y = (DATA_TYPE*)malloc(sizeof(DATA_TYPE)*N);
    z = (DATA_TYPE*)malloc(sizeof(DATA_TYPE)*N);

    for (i=0; i<N; i++) {
        y[i] = i * 1.0;
        z[i] = i * 2.0;
    }

    x_seq = 0;
    x_par = 0;


    // SEQUENTIAL CASE

    gettimeofday(&time, NULL);
    tstart_cpu = (double)clock()/CLOCKS_PER_SEC;
    tstart_wall = (double)time.tv_sec + (double)time.tv_usec * .000001;

    for (i=0; i<N; i++) x_seq += y[i] * z[i];

    tend_cpu = (double)clock()/CLOCKS_PER_SEC;
    gettimeofday(&time, NULL);
    tend_wall = (double)time.tv_sec + (double)time.tv_usec * .000001;


    cputime_par = tend_cpu - tstart_cpu;
    walltime_par = tend_wall - tstart_wall;
    cputime_par /= nthreads; // take the average cpu time per thread
    printf("Parallel CPU time  : %f\n", cputime_par);
    printf("Parallel Walltime  : %f\n", walltime_par);
    printf("Parallel result    : %f\n", x_par);

    // SPEEDUP

    printf("Speedup (cputime)  : %f\n", cputime_seq / cputime_par);
    printf("Speedup (walltime) : %f\n", walltime_seq / walltime_par);

    return 0;
}
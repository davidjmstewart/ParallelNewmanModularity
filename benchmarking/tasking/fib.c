#include <stdio.h>
#include <omp.h>
#include <time.h>

// gcc -fopenmp -g ./fib.c -o fib -lm -Wall -Wpedantic -Waggressive-loop-optimizations 
int fib(int n)
{
  
  int i, j;
  if (n<2)
    return n;
  else
    {
        // In case the sequence gets too short, execute the serial version
  if ( n < 20 )
  {
     return(fib(n-1) + fib(n -2 ));
  } else {
       #pragma omp task shared(i) firstprivate(n)
       i=fib(n-1);

       #pragma omp task shared(j) firstprivate(n)
       j=fib(n-2);

       #pragma omp taskwait
       return i+j;
    }
  }

}

int serialFib(int n)
{
  int i, j;
  if (n<2)
    return n;
  else
    {
       i=fib(n-1);

       j=fib(n-2);
       return i+j;
    }
}


int main()
{
    struct timespec start, finish;
    double elapsed;

  int n = 43;




    omp_set_dynamic(0);
    omp_set_num_threads(16);
    int parFib = 0;
    int seqFib = 0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    parFib = fib(n);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf ("Parallel took %f: %d\n", elapsed, parFib);


    clock_gettime(CLOCK_MONOTONIC, &start);
    seqFib = serialFib(n);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf ("Serial took %f: %d\n", elapsed, seqFib);

clock_gettime(CLOCK_MONOTONIC, &start);
  #pragma omp parallel shared(n)
  {
    // #pragma omp single
    printf ("fib(%d) = %d\n", n, fib(n));
  }

      clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf ("Par 2 took %f\n", elapsed);

}

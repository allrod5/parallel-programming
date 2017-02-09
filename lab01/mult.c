#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char **argv) {
   
    struct timeval tv;
    gettimeofday(&tv, 0);
    long t1 = tv.tv_sec*1000 + tv.tv_usec/1000;

    long **a, **b, **c;
    int N = 1000;

    if (argc == 2) {
      N = atoi (argv[1]);
      assert (N > 0);
    }

    int i,j,k,mul=2;
    int col_sum = N * (N-1) / 2;

    a = (long **)malloc (N * sizeof(long *));
    b = (long **)malloc (N * sizeof(long *));
    c = (long **)malloc (N * sizeof(long *));
    for (i=0; i<N; i++) {
      a[i] = (long *)malloc (N * sizeof(long));
      b[i] = (long *)malloc (N * sizeof(long));
      c[i] = (long *)malloc (N * sizeof(long));
    }


    for (i=0; i<N; i++)
      for (j=0; j<N; j++) {
	a[i][j] = i*mul;
	b[i][j] = i;
	c[i][j] = 0;
      }

    gettimeofday(&tv, 0);
    long t2 = tv.tv_sec*1000 + tv.tv_usec/1000;
    printf ("Matrix generation finished %ldms.\n", t2-t1);         
    
    #pragma omp parallel for private(j,k)
    for (i=0; i<N; i++)
      for (j=0; j<N; j++)
	for (k=0; k<N; k++)
	  c[i][j] += a[i][k] * b[k][j];

    gettimeofday(&tv, 0);
    long t3 = tv.tv_sec*1000 + tv.tv_usec/1000;
    printf ("Multiplication finished. %ldms.\n", t3-t2);
  
    for (i=0; i<N; i++)
      for (j=0; j<N; j++)
	assert ( c[i][j] == i*mul * col_sum);  

    gettimeofday(&tv, 0);
    long t4 = tv.tv_sec*1000 + tv.tv_usec/1000;
    printf ("Test finished. %ldms.\n", t4-t3);
}

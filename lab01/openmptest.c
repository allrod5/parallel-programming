#include <stdio.h>
#include <omp.h>

int main(int argc, char **argv) {

    omp_set_num_threads(4);
    int i, soma = 0;
    int a[] = {1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16};
    int b[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
    int c[] = {0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0};

    #pragma omp parallel
    {
        printf("1:Executing thread %d.\n", omp_get_thread_num());

        #pragma omp for
        for (i = 0; i < 16; i++) {
            c[i] += a[i] + b[i];
            printf( "t:%d i:%2d.\n", omp_get_thread_num(), i);
        }

        printf("1:Executing thread %d.\n", omp_get_thread_num());
    }

    for (i = 0; i < 16; i++)
        printf( "%d ", c[i]);
    
    printf("\n");

    return 0;
}

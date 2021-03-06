#include <stdio.h>

int main(int argc, char **argv) {
	omp_set_num_threads(4);
	int i, soma = 0;
	int a[] = {1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16};
	int b[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
	int c[] = {0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0};

	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (i = 0; i < 16; i++) {
			c[i] = a[i] + b[i];
			printf( "t:%d i:%2d.\n", omp_get_thread_num(), i);
		}
	}

	for (i = 0; i < 16; i++)
                printf( "%d ", c[i]);

        printf("\n\n");


        #pragma omp parallel
        {
                #pragma omp for schedule(dynamic)
                for (i = 0; i < 16; i++) {
                        c[i] = a[i] + b[i];
                        printf( "t:%d i:%2d.\n", omp_get_thread_num(), i);
                }
        }
	
	for (i = 0; i < 16; i++)
		printf( "%d ", c[i]);
	
	printf("\n");
}

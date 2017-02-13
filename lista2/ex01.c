#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

int main(int argc, char **argv) {
    int N, i, j, k;
    double desvio, grandeMedia, *media, **matriz;

    printf("Enter with the number of individuals: ");
    scanf("%d", &N);

    media = (double *) malloc(N * sizeof(double));

    printf("Enter with the number of tests per individual: ");
    scanf("%d", &k);

    matriz = (double **) malloc(N * sizeof(double *));
    for (i=0; i<N; i++) {
        matriz[i] = (double *) malloc(k * sizeof(double));
        for (j=0; j<k; j++) {
            printf("Enter the result of test %d for inividual %d: ", j, i);
            scanf("%lf", &matriz[i][j]);
        }
    }

    grandeMedia = 0;
    desvio = 0;

    #pragma omp parallel
    {
        #pragma omp for private(j) reduction(+:media)
        for (i = 0; i < N; i++)
            for (j = 0; j < k; j++)
                media[i] += matriz[i][j];

        #pragma omp barrier

        #pragma omp for reduction(/:media) reduction(+:grandeMedia)
        for (i = 0; i < N; i++) {
            media[i] /= k;
            grandeMedia += media[i];
        }
        grandeMedia /= N;

        #pragma omp barrier

        #pragma omp for reduction(+:desvio)
        for (i = 1; i < N; i++)
            desvio += (media[i] - grandeMedia) * (media[i] - grandeMedia);
    }

    desvio = sqrt(desvio/N);

    for (i=0; i<N; i++) {
        printf("\nThe mean score of individual %d was %lf", i, media[i]);
    }

    printf("\n\nThe mean score of all individuals was %lf\n", grandeMedia);
    printf("The standard deviation was %lf", desvio);

    return 0;
}
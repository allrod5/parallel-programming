#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <memory.h>

#define NUM_THREADS 4

typedef struct {
    int x;
    int y;
} int2;

typedef struct {
    double x;
    double y;
} double2;

typedef struct {
    int *i;
    int *progress;
    float **mat;
    double dt;
    int size;
    int ite;
    int nColumns;
    int nLines;
    pthread_mutex_t *i_mutex;
    pthread_mutex_t *progress_mutex;
} args;

int2 coordinatesConversion(double x,double y, int nColumns,int nLines) {

    int2 ret;

    //error return code=================================================
    int2 retError;
    retError.x=-1;
    retError.y=-1;
    //end===============================================================



    ret.x= (int) round(((2.0 + x) / 3.5) * ((double)(nColumns - 1)));
    ret.y= (int) round(((1.5 + y) / 3.5) * ((double)(nLines - 1)));

    //invalid parameters for  x or y arguments==========================
    if(ret.x<0 || ret.x>=nColumns) return retError;
    if(ret.y<0 || ret.y>=nLines) return retError;
    //end===============================================================
    return ret;
}

int printMatrixToFilePGM(float **mat, int tamx, int nLines, char *srcFile) {
    
    FILE *arq=fopen(srcFile,"w");
    
    int cont, cont2;
    float min,max;
    min=mat[0][0];
    max=mat[0][0];
    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<tamx;cont2++){
            if(min>mat[cont][cont2]) min=mat[cont][cont2];
            if(max<mat[cont][cont2]) max=mat[cont][cont2];
        }
    }
    max= (float) (max * 0.35);
    float delta=max-min;
    fprintf(arq,"P2 \n");
    fprintf(arq,"#comentario qualquer \n");
    fprintf(arq,"%d\n%d \n",tamx,nLines);
    fprintf(arq,"255\n");
    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<tamx;cont2++){
            int valpixel= (int) (((mat[cont][cont2] - min) / delta) * 255.0f);
            if(valpixel>255) valpixel=255;
            fprintf(arq,"%d \n", valpixel);
        }
    }
    fclose(arq);
    
    return 0;
}

float** mallocFloatMatrix(int tamx, int nLines, float defaultValueOfTheElementsAtMatrix) {
    float **errorCodeReturn=0x0;
    float **mat;
    float *data;
    int i, j;
    int condErrorMalloc=0; // error indicator at malloc
    
    // Contiguous multiarray
    data = malloc(sizeof (float) * tamx * nLines);
    mat=malloc(sizeof(float *)*nLines);
    if(mat==0x0) {
        printf("mallocFloatMatrix: malloc error: mat==0x0");
        return errorCodeReturn; // error at malloc return null pointer
    }
    for (i=0; i<nLines; i++) {
        mat[i] = &(data[tamx*i]);
    }


    for (i=0; i<nLines; i++) { // detect error at malloc
        if (mat[i]==0x0) {
            printf("mallocFloatMatrix: malloc error: mat[%i]==0x0", i);
            condErrorMalloc=1;
            break;
        }
    }

    if(condErrorMalloc==0){
        return mat;
    }
    for(i=0;i<nLines;i++){
        for(j=0;j<tamx;j++)
            mat[i][j]=defaultValueOfTheElementsAtMatrix;
    }
    for(i=0;i<nLines;i++)
        if(mat[i]!=0x0) free(mat[i]);

    free(mat);

    return errorCodeReturn;
}

int iteration(double x,double y, int nColumns,int nLines, int ite,int2 *iterationPath) {
    int cont;
    int condInvalidPointer=1;
    double2 z;
    z.x=0.0;
    z.y=0.0;
    double2 c;
    c.x=x;
    c.y=y;
    double2 zt;

    for(cont=0;cont<ite;cont++){
        //This  does z=z^2+c==================================================
        zt.x=((z.x*z.x)-(z.y*z.y))+c.x;
        zt.y=(2.0*(z.x*z.y))+c.y;
        z=zt;
        //end=================================================================
        if(((z.x*z.x)+(z.y*z.y))>4.0){
            if(cont>100)
                condInvalidPointer=0;
            break;
        }
        iterationPath[cont]=coordinatesConversion(z.x,z.y,nColumns,nLines);
    }
    if(condInvalidPointer)
        return 0;

    return cont;
}

int main(int argc, char *argv[]) {
    //image size=========================================================
    int nColumns = 4096;
    int nLines = 4096;
    //end================================================================
    //size points========================================================
    double dt = 0.00025; // quantity of points going to increase with the decrease of the dt value
    int size = (int) round(4.0 / dt); // sizeOfPoints = size*size
    //end================================================================
    int ite = 600;
    
    float **mat = mallocFloatMatrix(nColumns, nLines, 0.0f);
    if (mat == 0x0) {
        return 0;
    }
    
    int i, j, k;
    double x, y;
    int progress=0;

    int numprocs, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int start = rank * (size/numprocs + (size%numprocs-rank > 0 ? 1 : 0));
    int end = (rank+1) * (size/numprocs + (size%numprocs-(rank+1) > 0 ? 1 : 0));

    omp_set_num_threads(5);
    
    #pragma omp parallel for private(k, x, y)
    for (i=start; i<end; i++) { // real component of C at $z_{n+1}=z_n+C$
        
        x = -2.0 + i*dt;
        for (y=-2.0; y<2.0; y=y+dt) { // imaginary component of C at $z_{n+1}=z_n+C$

            int2* iterationPath = (int2 *) malloc(sizeof(int2)*ite);
            if (iterationPath==0x0) {
                //return 0x0;
                i=end;
            }

            // completedIterations = quantity of elements at vector iterationPath
            int completedIterations = iteration(x, y, nColumns, nLines, ite, iterationPath);
            for (k=0; k<completedIterations; k++) {
                if (iterationPath[k].x!=-1 && iterationPath[k].y!=-1) { // test if a point z in the iteration k may be normalized to coordinates at matrix mat.
                    mat[iterationPath[k].x][iterationPath[k].y] = mat[iterationPath[k].x][iterationPath[k].y] + 1.0f; // increments a point in matrix, this point is pointed by z with  z points normalized.
                }
            }
            free(iterationPath);
        }

        progress++;
        if (progress%100 == 0) { // print at screen information about progress of the operation
            printf("%lf \n", x);
        }
    }
    
    if (rank!=0) {
        MPI_Send(&(mat[0][0]), nColumns * nLines, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
    } else {
        float **aux;
        for (i=1; i<numprocs; i++) {
            aux = mallocFloatMatrix(nColumns, nLines, 0.0f);
            MPI_Status status;
            MPI_Recv(&(aux[0][0]), nColumns * nLines, MPI_FLOAT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            for (j = 0; j < nLines; j++) {
                for (k = 0; k < nColumns; k++) {
                    mat[j][k] += aux[j][k];
                }
            }
            free(aux[0]);
            free(aux);
        }
        
        printMatrixToFilePGM(mat, nColumns, nLines, "output.pgm");
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    free(mat[0]);
    free(mat);
    
    return 0;
}
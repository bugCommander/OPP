#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N  1000
#define t 10e-5
#define e 10e-8
///#define _ALLOCATION_LIMIT 4
int __allocs = 0;

static void *my_malloc(size_t size) {
#ifdef _ALLOCATION_LIMIT
    if (__allocs > _ALLOCATION_LIMIT) return NULL;
#endif
    __allocs++;
    return malloc(size);
}

static void my_free(void *ptr) {
    __allocs--;
    free(ptr);
}

#define malloc my_malloc
#define free my_free
void my_exit(int err) {
    fprintf(stderr, "Allocation balance %d %s\n", __allocs, __allocs ? "- Memory leaks occur!!!" : "is ok");
    exit(err);
}
#define exit my_exit


double abs_vector(const double *vector){
    double vector_size = 0;
    for(size_t i = 0;i < N;++i){
        vector_size += vector[i]*vector[i];
    }
    return sqrt(vector_size);
}


void matrix_vector_mult( double * matrix,  double * vector, double*res){

    for(size_t i = 0;i < N;++i) {
        for (size_t j = 0; j < N; ++j) {
            res[i] += matrix[i + N * j] * vector[j];
        }
    }

}


void scalar_mult_vector( double * vector, double scalar, double * res){
    for(size_t i = 0; i < N;++i){
        res[i] = vector[i] * scalar;
    }

}
void vector_sub(const double *A, const double *B,double *res){
    for(size_t i = 0;i < N;++i){
        res[i] = A[i] - B[i];
    }
}
int main(int argc, char **argv) {
    ///MPI_Init(&argc,&argv);
    double *A = (double *)malloc(sizeof(double )*N*N);
    double *b = (double *)malloc(sizeof(double)*N);
    double *res = (double*) malloc(sizeof(double) * N);
    for(int i = 0;i < N;++i){
        b[i] = N+1;
        res[i] = 0;
        for(int j = 0; j < N;++j){
           if(i == j){
               A[i + N*j] = 2.0;
           }
           else {
               A[i + N*j] = 1.0;

           }
        }
    }
    double *aux_vector = (double*) malloc(sizeof(double) * N); ///Ax^n - b
    while(1){
        matrix_vector_mult(A,res,aux_vector); ///Ax^n
        vector_sub(aux_vector,b,aux_vector); ///Ax^n - b
        if(abs_vector(aux_vector) < e){
            break;
        }
        scalar_mult_vector(aux_vector,t,aux_vector);
        vector_sub(res,aux_vector,res);



    }
    for(int i = 0;i < N;++i){
            printf("%f \t",res[i]);
        }
        printf("\n\t");

    free(b);
    free(A);
    free(res);
    free(aux_vector);
    exit(0);
}
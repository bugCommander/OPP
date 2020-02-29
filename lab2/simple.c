#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N  1000
#define t 10e-5
#define e 10e-8
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
    double *A = (double *)malloc(sizeof(double )*N*N);
    double *b = (double *)malloc(sizeof(double)*N);
    double *x = (double*) malloc(sizeof(double) * N);
    for(int i = 0;i < N;++i){
        b[i] = N+1;
        x[i] = 0;
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
        matrix_vector_mult(A, x, aux_vector); ///Ax^n
        vector_sub(aux_vector,b,aux_vector); ///Ax^n - b
        if(abs_vector(aux_vector)/abs_vector(b) < e){
            break;
        }
        scalar_mult_vector(aux_vector,t,aux_vector);
        vector_sub(x, aux_vector, x);



    }
    for(int i = 0;i < N;++i){
        printf("%f \t", x[i]);
    }
    printf("\n\t");

    free(b);
    free(A);
    free(x);
    free(aux_vector);
    exit(0);
}


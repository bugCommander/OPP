//
// Created by ksandr on 06.03.2020.
//
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N  100
#define t 10e-6
#define e 10e-9
#define type guided ////static  or dynamic
int main(int argc, char **argv) {
/*
    if(argc != 2){
        printf("wrong param");
        exit(0);
    }
    */

     omp_set_num_threads(2);
   /// omp_set_num_threads(atoi(argv[1]));

    double *matrix = (double *) malloc(sizeof(double) * N * N);
    double *b = (double *) malloc(sizeof(double) * N);
    double *x = (double *) malloc(sizeof(double) * N);
    double *tmp_x = (double *) malloc(sizeof(double) * N);
    int i, j;

    for (i = 0; i < N; ++i) {
        b[i] = N + 1;
        x[i] = 0;
        for (j = 0; j < N; ++j) {
            if (i == j) {
                matrix[i * N + j] = 2.0;
            } else {
                matrix[i * N + j] = 1.0;

            }
        }
    }
    double b_length = 0;
    double t1 = omp_get_wtime();
//#pragma omp parallel for reduction(+:b_length)
#pragma omp for schedule  (type)
    for (i = 0; i < N; ++i) {
        b_length += b[i] * b[i];
    }
    b_length = sqrt(b_length);

    while (1) {
        double result_length = 0;

#pragma omp parallel for schedule(type)
        for (int i = 0; i < N; ++i) {
            double aux_sum = 0;
#pragma omp parallel for schedule(type)
            for (int j = 0; j < N; ++j) {
                aux_sum += matrix[i * N + j] * x[j]; //Ax
            }
            aux_sum = aux_sum - b[i]; //Ax - b
            tmp_x[i] = x[i] - aux_sum * t; //x^n=x-t(Ax-b)

            result_length += aux_sum * aux_sum;
        }
#pragma omp  single
        for (int i = 0; i < N; ++i) {
            x[i] = tmp_x[i];
        }


        result_length = sqrt(result_length);
        if(result_length/b_length <e){
            break;
        }

    }
    double t2 = omp_get_wtime() - t1;

    for (int i = 0; i < N; ++i) {
        printf("%f \t", x[i]);
    }

    printf("\n\t");
    printf("%f\n", t2);
    free(tmp_x);
    free(b);
    free(matrix);
    free(x);
    exit(0);
}


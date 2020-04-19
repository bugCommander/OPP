//
// Created by ksandr on 17.04.2020.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openmpi/mpi.h>

#define N 100
#define a 1e5

#define e 1e-8
const double min_x = -1;
const double min_y = -1;
const double min_z = -1;
const double max_x = 1;
const double max_y = 1;
const double max_z = 1;

double phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double ro(double x, double y, double z) {
    return 6 - a * phi(x, y, z);
}

void init_matrices(double *matrix, double *aux_matrix,
                   int start_i, int proc_size,
                   double delta_x, double delta_y, double delta_z) {
    double cur_x, cur_y, cur_z;
    int real_i;
    int pos;

    for (int i = 0; i < proc_size; ++i) {
        real_i = start_i + i;
        cur_z = min_z + real_i * delta_z;


        for (int j = 0; j < N; ++j) {
            cur_x = min_x + j * delta_x;
            for (int k = 0; k < N; ++k) {
                cur_y = min_y + k * delta_y;
                pos = i * N * N + j * N + k;
                if (real_i == 0 || real_i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                    matrix[pos] = phi(cur_x, cur_y, cur_z);
                } else {
                    matrix[pos] = 0;
                }
                aux_matrix[pos] = matrix[pos];

            }
        }

    }


}

double max(double x, double y) {
    return x >= y ? x : y;
}

double recalc_layer(double *matrix, double *aux_matrix,
                    int start_i, int i,
                    double delta_x, double delta_y, double delta_z,
                    double constant) {
    double real_i = start_i + i;
    if (real_i == 0 || real_i == N - 1) {
        return 0;
    }
    double max_delta = 0;
    double real_z = min_z + delta_z * real_i;
    double real_x, real_y;
    int pos;
    for (int j = 0; j < N; ++j) {
        real_x = min_x + delta_x * j;
        for (int k = 0; k < N; ++k) {
            if (j==0 || j==N-1 || k==0 || k==N-1) {
                aux_matrix[i*N*N+i*N+j] = matrix[i*N*N+i*N+j];
                continue;
            }
            real_y = min_y + delta_y * k;
            double x = (aux_matrix[i * N * N + (j - 1) * N + k] + aux_matrix[i * N * N + (j + 1) * N + k]) /
                       (delta_x * delta_x);
            double y = (aux_matrix[i * N * N + j * N + (k - 1)] + aux_matrix[i * N * N + j * N + (k + 1)]) /
                       (delta_y * delta_y);
            double z = (aux_matrix[(i - 1) * N * N + j * N + k] + aux_matrix[(i + 1) * N * N + j * N + k]) /
                       (delta_z * delta_z);
            pos = i * N * N + j * N + k;
            matrix[pos] = constant * (x + y + z - ro(real_x, real_y, real_z));
            max_delta = max(fabs(matrix[pos] - aux_matrix[pos]), max_delta);
        }
    }
    return max_delta;


}

int main(int argc, char **argv) {

    double delta_x = (max_x - min_x) / (N - 1);
    double delta_y = (max_y - min_y) / (N - 1);
    double delta_z = (max_z - min_z) / (N - 1);
    double constant = 1 / (2 / (delta_x * delta_x) + 2 / (delta_y * delta_y) + 2 / (delta_z * delta_z) + a);
    int size, rank, proc_size;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (N % size) {
        perror("can't create the world");
        MPI_Finalize();
        return -1;
    }
    proc_size = (rank == size - 1 || rank == 0) ? N / size + 1 : N / size + 2;
    int start_i = rank*(N/size) -1;



    double *matrix_part = (double *) malloc(sizeof(double) * N * N * (proc_size));
    double *aux_matrix_part = (double *) malloc(sizeof(double) * N * N * proc_size);
    init_matrices(matrix_part, aux_matrix_part, start_i, proc_size, delta_x, delta_y, delta_z);

    MPI_Request request_arr[4];
    while (1) {
        double max_delta;
        double tmp_max_delta;
        max_delta = recalc_layer(matrix_part, aux_matrix_part, start_i, 1, delta_x, delta_y,
                                 delta_z, constant);/// upper layer
        tmp_max_delta = recalc_layer(matrix_part, aux_matrix_part, start_i, proc_size - 2, delta_x,
                                     delta_y, delta_z, constant); /// lower layer
        max_delta = max(tmp_max_delta, max_delta);
        if (rank != 0) {
            MPI_Isend(matrix_part + N * N, N * N, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &request_arr[0]);
            MPI_Irecv(matrix_part, N * N, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &request_arr[2]);
        }
        if (rank != size - 1) {
            MPI_Isend(matrix_part + (proc_size - 2) * N * N, N * N, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD,
                      &request_arr[1]);
            MPI_Irecv(matrix_part + (proc_size - 1) * N * N, N * N, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD,
                      &request_arr[3]);

        }

        if (rank != 0) {
            MPI_Wait(&request_arr[0], MPI_STATUSES_IGNORE);
            MPI_Wait(&request_arr[2], MPI_STATUS_IGNORE);

        }

        if (rank != size - 1) {
            MPI_Wait(&request_arr[1], MPI_STATUS_IGNORE);
            MPI_Wait(&request_arr[3], MPI_STATUS_IGNORE);

        }

        for (int i = 2; i < proc_size - 2; ++i) {
            tmp_max_delta = recalc_layer(matrix_part, aux_matrix_part, start_i, i, delta_x,
                                         delta_y, delta_z, constant);
            max_delta = max(tmp_max_delta, max_delta);
        }
        for (int i = 0; i < proc_size; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 1; k < N; ++k) {
                    int pos = i * N * N + j * N + k;
                    aux_matrix_part[pos] = matrix_part[pos];
                }
            }
        }
        double reduce_max_delta = 0;
        MPI_Reduce(&max_delta, &reduce_max_delta, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&reduce_max_delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (reduce_max_delta < e) {
            break;

        }
    }

/*
    printf("%d\n", rank);
    for (int i = 0; i < proc_size; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                printf("\tmatrix[%d][%d][%d] == %lf \t", j, k, start_i + i, matrix_part[i * N * N + N * j + k]);
            }
            printf("\n");
        }
        printf("\nz_step()\n");

    }
    */


    free(matrix_part);
    free(aux_matrix_part);
    MPI_Finalize();
    return 0;
}
//
// Created by ksandr on 14.02.2020.
//
#include <stdio.h>
#include <openmpi/mpi.h>
#include <stdlib.h>
#include <math.h>

#define N 4
#define t 10e-6
#define e 10e-9
#define length_tag 42
#define x_tag 43
#define b_length_tag 44

int main(int argc, char **argv) {


    int rank, world_rank;
    double start;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_rank);
    int multiple_processes = world_rank - (N % world_rank); /// число процессов в которых кратное N  число строк

    int *lines_counter = (int *) malloc(sizeof(int) * world_rank);
    for (int i = 0; i < world_rank; ++i) {
        if (i < multiple_processes) {
            lines_counter[i] = N / world_rank;
        } else {
            lines_counter[i] = N / world_rank + 1;
        }

    }


    double *matrix = (double *) malloc(N * lines_counter[rank] * sizeof(double));
    int skip_lines = 0;
    for (int i = 0; i < rank; ++i) {
        skip_lines += lines_counter[i];
    }
    for (int i = 0; i < lines_counter[rank]; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i + skip_lines == j) {
                matrix[i * N + j] = 2;

            } else {
                matrix[i * N + j] = 1;

            }

        }
    }


    double *part_x = (double *) malloc(sizeof(double) * (lines_counter[rank]));
    double *part_b = (double *) malloc(sizeof(double) * lines_counter[rank]);
    double part_b_length = 0;
    for (int i = 0; i < lines_counter[rank]; ++i) {
        part_x[i] = 0;
        part_b[i] = N + 1;
        part_b_length += (part_b[i] * part_b[i]);
    }


    double *x = (double *) malloc(sizeof(double) * N);

    for (int i = 0; i < world_rank; ++i) {
        if (i != rank) {
            MPI_Send(part_x, lines_counter[rank], MPI_DOUBLE, i, x_tag, MPI_COMM_WORLD);
        }

    }


    short is_end = 0;
    if (rank == 0) {
        start = MPI_Wtime();
    }



    while (!is_end) {
        int current_line = 0;
        for(int i = 0; i < lines_counter[rank];++i){
            x[skip_lines+i] = part_x[i];
        }
        for (int i = 0; i < world_rank; ++i) {
            if (i != rank) {

                MPI_Recv(&x[current_line], lines_counter[i], MPI_DOUBLE, i, x_tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            }
            current_line += lines_counter[i];

        }

        double part_of_length = 0;
        for (int i = 0; i < lines_counter[rank]; ++i) {
            double aux_sum = 0; //// string_matrix * vector = 1 number ( this is)
            for (int j = 0; j < N; ++j) {
                aux_sum += matrix[i * N + j] * x[j];
            }
            aux_sum -= part_b[i];

            part_x[i] = x[skip_lines + i] - t * (aux_sum);
            part_of_length += aux_sum * aux_sum; ///  length of Ax^n - b


        }

        if (rank != 0) {
            MPI_Send(&part_of_length, 1, MPI_DOUBLE, 0, length_tag, MPI_COMM_WORLD);
            MPI_Send(&part_b_length, 1, MPI_DOUBLE, 0, b_length_tag, MPI_COMM_WORLD);
        }


        for (int i = 0; i < world_rank; ++i) {
            if (i != rank) {
                MPI_Send(part_x, lines_counter[rank], MPI_DOUBLE, i, x_tag, MPI_COMM_WORLD);
            }

        }

        if (rank == 0) {
            double length = part_of_length;
            double length_b = part_b_length;
            for (int i = 1; i < world_rank; ++i) {
                double r_length = 0;
                MPI_Recv(&r_length, 1, MPI_DOUBLE, i, length_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                length += r_length;
                MPI_Recv(&r_length, 1, MPI_DOUBLE, i, b_length_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                length_b += r_length;
            }


            if ((sqrt(length) / sqrt(length_b)) < e) {
                is_end = 1;
            }


        }
        MPI_Bcast(&is_end, 1, MPI_SHORT, 0, MPI_COMM_WORLD);


    }

    if (rank == 0) {

        double time = MPI_Wtime() - start;

        for (int j = 0; j < N; ++j) {
            printf("%f ", x[j]);
        }
        printf("\n");

        printf("time %f in %d processes", time, world_rank);
    }


    free(lines_counter);
    free(part_x);
    free(matrix);
    free(x);
    free(part_b);
    MPI_Finalize();


    return 0;

}



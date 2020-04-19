#include <stdio.h>
#include <stdlib.h>
#include <openmpi/mpi.h>
#include <string.h>
#include <math.h>

double max(double x, double y) {
    return x >= y ? x : y;
}

const double eps = 1e-8;
const double a = 1e5;
const int N = 100;
const double min_x = -1, min_y = -1, min_z = -1;
const double max_x = 1, max_y = 1, max_z = 1;

double phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double ro(double x, double y, double z) {
    return 6 - a * phi(x, y, z);
}

double recalc_layer(double *matrix, double *aux_matrix,
             int start_i, int i,
             double delta_x, double delta_y, double delta_z,
             double constant) {
    int real_z = start_i + i;
    if (real_z == 0 || real_z == N - 1) {///border OZ
        memcpy(aux_matrix + i * N * N, matrix + i * N * N, N * N * sizeof(double));
        return 0;
    }
    double max_delta = 0;
    double cur_z = min_z + real_z * delta_z;
    for (int j = 0; j < N; j++) {
        double cur_x = min_x + j * delta_x;
        for (int k = 0; k < N; k++) {
            double cur_y = min_y + k * delta_y;
            if (j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                aux_matrix[i * N * N + j * N + k] = matrix[i * N * N + j * N + k];
                continue;///border OX OY
            }
            int pos = i * N * N + j * N + k;
            //// asum  phi^(m+1)
            double x = (matrix[i * N * N + (j + 1) * N + k] + matrix[i * N * N + (j - 1) * N + k]) /
                       (delta_x * delta_x);
            double y = (matrix[i * N * N + j * N + (k + 1)] + matrix[i * N * N + j * N + (k - 1)]) /
                       (delta_y * delta_y);
            double z = (matrix[(i + 1) * N * N + j * N + k] + matrix[(i - 1) * N * N + j * N + k]) /
                       (delta_z * delta_z);
            aux_matrix[pos] = constant * (x + y + z - ro(cur_x, cur_y, cur_z));
            max_delta = max(max_delta, fabs(aux_matrix[pos] - matrix[pos]));
        }
    }
    return max_delta;
}

void fill_part_matrix(double *matrix_part, int proc_size, int start_i, double delta_x, double delta_y, double delta_z) {
    for (int i = 0; i < proc_size; i++) {
        int real_i = i + start_i;
        double cur_z = min_z + delta_z * real_i;
        for (int j = 0; j < N; j++) {
            double cur_x = min_x + delta_x * j;
            for (int k = 0; k < N; k++) {
                double cur_y = min_y + delta_y * k;
                int pos = i * N * N + j * N + k;
                if (real_i == 0 || real_i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                    matrix_part[pos] = phi(cur_x, cur_y, cur_z);
                } else {
                    matrix_part[pos] = 0;
                }
            }
        }
    }
}

void Icommutation(double *matrix_part, int proc_size, int rank, int size) {
    MPI_Request request_arr[4];
    if (rank != 0) {
        ///senging upper layer to prev proc( becoming lower layer)
        MPI_Isend(matrix_part + N * N, N * N, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &request_arr[0]);
        ///receiving lower layer ( becoming upper layer)
        MPI_Irecv(matrix_part, N * N, MPI_DOUBLE, rank - 1, 42, MPI_COMM_WORLD, &request_arr[1]);
    }
    if (rank != size - 1) {
        //// senging lower layer to becoming  upper layer
        MPI_Isend(matrix_part + proc_size * N * N, N * N, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD,
                  &request_arr[2]);
        ////receiving upper layer ( becoming lower layer)
        MPI_Irecv(matrix_part + (proc_size +1) * N * N, N * N, MPI_DOUBLE, rank + 1, 42, MPI_COMM_WORLD,
                  &request_arr[3]);

    }
}
double find_max_delta(double *matrix,int proc_size,int start_i, double delta_x, double delta_y, double delta_z){
    double max = 0;
    double real_x,real_y,real_z;
    for(int i = 0; i < proc_size;++i){
        real_z = min_z + (i+start_i)*delta_z;
        for(int j = 0; j < N;++j){
            real_x = min_x + j*delta_x;

            for(int k = 0; k < N; ++k){
                real_y = min_y + k*delta_y;

                max = (fabs(matrix[i*N*N+j*N+k] - phi(real_x,real_y,real_z)),max);

            }
        }
    }
    return max;

}
int main(int argc, char **argv) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double delta_x = (max_x - min_x) / (N - 1);
    double delta_y = (max_y - min_y) / (N - 1);
    double delta_z = (max_z - min_z) / (N - 1);
    double constant = 1 / (2 / (delta_x * delta_x) + 2 / (delta_y * delta_y) + 2 / (delta_z * delta_z) + a);
    int proc_size = N / size;
    int alloc_size = proc_size+2; //// for new local borders
    if (N % size) {
        if (rank == 0) {
            printf("N mod size should be 0");
        }
        MPI_Finalize();
        return 0;
    }
    double *matrix_part = (double *) malloc(sizeof(double) * alloc_size * N * N);
    if(matrix_part ==NULL){
        printf("can't allocate matrix in %d rank",rank);
        MPI_Finalize();
        return -2;
    }
    double *aux_matrix_part = (double *) malloc(sizeof(double) * alloc_size * N * N);
    if(aux_matrix_part == NULL){
        printf("can't allocate aux_matrix_part in %d rank",rank);
        free(matrix_part);
        MPI_Finalize();
        return -2;

    }
    int start_i = rank * proc_size - 1;




    fill_part_matrix(matrix_part, alloc_size, start_i, delta_x, delta_y, delta_z);
    double t1 = MPI_Wtime();
    while (1) {
        double max_delta;
        double tmp_delta;
        ///upper layer;
        max_delta = recalc_layer(matrix_part, aux_matrix_part, start_i, 1, delta_x, delta_y, delta_z, constant);
        ///lower layer;
        tmp_delta = recalc_layer(matrix_part, aux_matrix_part, start_i, proc_size, delta_x, delta_y, delta_z,
                                 constant);
        max_delta = max(max_delta, tmp_delta);
        Icommutation(aux_matrix_part, proc_size, rank, size);
        //// layers between lower & upper layers
        for (int i = 2; i < proc_size; i++) {

            max_delta = max(max_delta,
                            recalc_layer(matrix_part, aux_matrix_part, start_i, i, delta_x, delta_y, delta_z,
                                         constant));
        }
        memcpy(matrix_part, aux_matrix_part, alloc_size * N * N * sizeof(double));
        double max_delta_shared;
        MPI_Reduce(&max_delta, &max_delta_shared, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_delta_shared, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (max_delta_shared < eps) {
            break;
        }
    }
    double max_delta = find_max_delta(matrix_part,proc_size,start_i,delta_x,delta_y,delta_z);
    double max_delta_shared;
    MPI_Reduce(&max_delta, &max_delta_shared, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    double *res;
    int flag = 0;
    if (rank == 0) {
        res = (double*) malloc(sizeof(double)*N*N*N);
        if(res == NULL){
            printf("can't allocate res");
            flag =1;
            MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);

        }else {
            printf("time: %lf\ndelta : %lf\n", MPI_Wtime() - t1, max_delta_shared);
        }
    }
    if(!flag) {
        MPI_Gather(matrix_part, N * N * proc_size, MPI_DOUBLE, res, N * N * proc_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //// builging  function field

    }
    if(!flag && rank == 0){
        free(res);

    }
    free(aux_matrix_part);
    free(matrix_part);
    MPI_Finalize();
    return 0;
}
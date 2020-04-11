#include <stdio.h>
#include <stdlib.h>
#include <openmpi/mpi.h>

#define A_rows 16 ///n1
#define A_cols 8 ///n2
#define B_rows 8 /// == n2
#define B_cols 16 /// == n3
#define grid_columns 2
#define grid_rows 2
#define n_dims 2
#define x 0
#define y 1

void init_grid_comm(MPI_Comm *grid_comm, int coords[n_dims], int dims[n_dims]) {
    int periods[n_dims] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, n_dims, dims, periods, 0, grid_comm);
    int rank;
    MPI_Comm_rank(*grid_comm, &rank);
    MPI_Cart_coords(*grid_comm, rank, n_dims, coords);
}


void init_helper_comms(MPI_Comm grid, MPI_Comm *rows_comm, MPI_Comm *col_comm, int coords[n_dims]) {
    MPI_Comm_split(grid, coords[y], coords[x], rows_comm);
    MPI_Comm_split(grid, coords[x], coords[y], col_comm);
}


void distribute_matrices(double *matrixA, double *matrixB, double *partA, double *partB, MPI_Comm row_comm,
                         MPI_Comm col_comm,const int coords[n_dims]) {
    if(coords[x] == 0) {
        MPI_Scatter(matrixA, A_rows * A_cols / grid_rows, MPI_DOUBLE, partA, A_rows * A_cols / grid_rows, MPI_DOUBLE, 0,
                    col_comm);
    }
    MPI_Bcast(partA, A_rows * A_cols / grid_rows, MPI_DOUBLE, 0, row_comm);
    if(coords[y] == 0){

        MPI_Datatype aux_t;
        MPI_Datatype b_column_t;
        MPI_Type_vector(B_rows, B_cols / grid_columns, B_cols, MPI_DOUBLE, &aux_t);
        MPI_Type_commit(&aux_t);
        MPI_Type_create_resized(aux_t, 0, B_cols/ grid_columns * sizeof(double), &b_column_t);
        MPI_Type_commit(&b_column_t);

        MPI_Datatype recv_t;
        PMPI_Type_contiguous(B_rows * B_cols / grid_columns, MPI_DOUBLE, &recv_t);
        MPI_Type_commit(&recv_t);

        MPI_Scatter(matrixB, 1, b_column_t,partB, 1, recv_t,0, row_comm);

        MPI_Type_free(&b_column_t);
        MPI_Type_free(&aux_t);
        MPI_Type_free(&recv_t);
    }

    MPI_Bcast(partB,B_rows * B_cols / grid_columns,MPI_DOUBLE,0,col_comm);



}

void fill_matrix(double *matrix, int row, int col, double value) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (i == j) {
                matrix[i * col + j] = 4;
            } else {
                matrix[i * col + j] = value;
            }
        }
    }
}


void matrix_mult(const double *A, const double *B, double *C,int ARows,int AColumns, int BRows,int BColumns){
    for (int i = 0; i < ARows; i++) {
        for (int j = 0; j < BColumns; j++) {
            double sum = 0;

            for (int k = 0; k < AColumns; k++) {
                sum += A[i * AColumns + k] * B[k * BColumns + j];
            }

            C[i * BColumns + j] = sum;
        }
    }
}




int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int dims[n_dims] = {grid_columns, grid_rows};
    int coords[n_dims];
    MPI_Comm grid_comm, row_comm, col_comm;
    init_grid_comm(&grid_comm, coords, dims);
    init_helper_comms(grid_comm, &row_comm, &col_comm, coords);
    int rank;
    MPI_Comm_rank(grid_comm, &rank);
    double *A = NULL;
    double *B = NULL;

    double *C =NULL;
    double *partC = (double*)malloc(sizeof(double)*A_rows*B_cols/(grid_columns*grid_rows));
    double *partA = (double *) malloc(sizeof(double) * A_rows * A_cols / grid_rows);
    double *partB = (double *) malloc(sizeof(double) * B_rows * B_cols / grid_columns);


    if (rank == 0) {
        A = (double *) malloc(sizeof(double) * A_rows * A_cols);
        B = (double*) malloc(sizeof(double)*B_rows*B_cols);
        C = (double*)malloc(sizeof(double)*A_rows*B_cols);
        fill_matrix(A, A_rows, A_cols, 2);
        fill_matrix(B,B_rows,B_cols,2);
        printf("matrix A\n");
        for (int i = 0; i < A_rows; ++i) {
            for (int j = 0; j < A_cols; ++j) {
                printf("%lf ", A[i * A_cols + j]);
            }
            printf("\n");
        }
        printf("matrix B\n");


        for (int i = 0; i < B_rows; ++i) {
            for (int j = 0; j < B_cols; ++j) {
                printf("%lf ", B[i * B_cols + j]);
            }
            printf("\n");
        }
        printf("//////////////\n");




    }
    distribute_matrices(A,B,partA,partB,row_comm,col_comm,coords);
    matrix_mult(partA,partB,partC,A_rows/grid_rows,A_cols,B_rows/grid_rows,B_cols/grid_columns);
    MPI_Gather(partC,A_rows*B_cols/(grid_columns*grid_rows),MPI_DOUBLE,C,A_rows*B_cols/(grid_columns*grid_rows),MPI_DOUBLE,0,grid_comm);
    if(rank ==0){
        for(int i =0; i < A_rows;++i) {
            for (int j = 0; j < B_rows; ++j) {
                printf("%lf ", C[i * B_rows + j]);
            }
            printf("\n");
        }



    }


/*
    printf(" coords is %d %d \n", coords[x], coords[y]);
    for(int i = 0; i < A_rows/grid_rows;++i){
        for(int j = 0; j < B_cols/grid_columns;++j){
            printf("%lf ",partC[i*B_cols/grid_columns+j]);
        }
        printf("\n");
    }
    /*
    /*
        for (int i = 0; i < B_rows / grid_rows; ++i) {
            for (int j = 0; j < A_cols; ++j) {
                printf("%lf ", partA[i * A_cols + j]);
            }
            printf("\n");
        }
        */
    /*
    for(int i = 0; i < B_rows;++i){
        for(int j = 0; j < B_cols/grid_columns;++j)
        printf("%lf ", partB[i*B_cols/grid_columns+j]);
    }
    printf("\n");
     */
    if(rank ==0){
        free(A);
        free(B);
        free(C);
    }
    free(partA);
    free(partB);
    free(partC);




    MPI_Comm_free(&grid_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();
    return 0;
}

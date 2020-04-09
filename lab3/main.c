#include <stdio.h>
#include <openmpi/mpi.h>
#include <stdlib.h>

#define nrows 2
#define ncols 2
#define ndims 2
#define n1 2
#define n2 4
#define n3 4

double *fillMatrix(int row, int col, double value) {
    double *matrix = (double *) malloc(sizeof(double) * row * col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
                if(i!=j)
                matrix[i * col + j] = value;
                else matrix[i * col + j] = 4;
            }

        }



    return matrix;
}


int main(int argc, char **argv) {
    int dims[ndims] = {nrows, ncols};
    int periods[ndims] = {0, 0};
    int coords[ndims];
    int reorder = 0;
    ///double *partA, partB, PartC;
    int A_rows = n1/nrows;
    int B_cols = n3/ncols;



    int rank;
    double *A, *B, *C;
    MPI_Comm gridComm;
    MPI_Comm rowComm, collComm;

    MPI_Init(&argc, &argv);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &gridComm);
    MPI_Comm_rank(gridComm, &rank);
    MPI_Cart_coords(gridComm, rank, ndims, coords);
    if (rank == 0) {
        A = fillMatrix(n1,n2,1.0);
        B = fillMatrix(n1,n2,2.0);
       C = fillMatrix(n1,n3,0);
        /*
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                printf("%lf ", A[i * n2 + j]);

            }
            printf("\n");


        }
         */


    }
    ///color - coords[0] - row - unique color /// key - coords [1] - ordering
    MPI_Comm_split(gridComm, coords[0], coords[1], &rowComm);
    MPI_Comm_split(gridComm, coords[1], coords[0], &collComm);
    double *partA = (double*)malloc(sizeof(double)*A_rows*n2);
    double *partB = (double*)malloc(sizeof(double)*n2*B_cols);
    double *partC = (double*)malloc(sizeof(double)*A_rows*B_cols);

    /// тут надо раскидать по столбцам и строкам матрицы
    ///
    if (coords[0] == 0) {
        MPI_Scatter(B,B_cols*n2,MPI_DOUBLE,partB,B_cols*n2,MPI_DOUBLE,0,rowComm);



    }
    MPI_Bcast(partB,n2*B_cols,MPI_DOUBLE,0,collComm);

    if (coords[1] == 0) {
        MPI_Scatter(A,A_rows*n2,MPI_DOUBLE,partA,A_rows*n2,MPI_DOUBLE,0,collComm);
        ////send &recv a in the row




    }
    MPI_Bcast(partA,n2*A_rows,MPI_DOUBLE,0,rowComm);
    if (rank == 0) {


    }
    printf("coordinates is %d %d\n",coords[0],coords[1]);
    for(int i = 0;i < n3;++i){
        for(int j = 0; j < B_cols;++j){
            printf("%lf\n",partB[i*B_cols+j]);

        }
        printf("\n");
    }


    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&collComm);
    MPI_Comm_free(&gridComm);
    MPI_Finalize();


    return 0;
}

#include <stdio.h>
#include <openmpi/mpi.h>
#include <stdlib.h>

#define n1 2
#define n2 2
#define n3 2
void read_matrix(double *matrix,int row,int colon,FILE*input){
    for(int i = 0;i < row;++i){
        for(int j =0; j < colon; ++j)
            fscanf(input,"%lf",&matrix[i*row + j]);
    }
    
}


void write_matrix(double *matrix,int row,int colon,FILE*output){
    for(int i = 0;i < row;++i){
        for(int j =0; j < colon; ++j) {
            fprintf(output,"%lf ", matrix[i * row + j]);
        }
        fprintf(output,"\n");
    }
}
int main(int argc, char **argv) {
    int world_size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) { //// считывание матрицы
        FILE *input = fopen("A.txt", "r");
        if (!input) {
            MPI_Finalize();
            return -1;
        }
        double *matrix_A = (double *) malloc(sizeof(double) * n1 * n2);
        double *matrix_B = (double *) malloc(sizeof(double) * n2 * n3);
        double *matrix_C = (double *) malloc(sizeof(double) * n1 * n3);
        read_matrix(matrix_A,n1,n2,input);
        fclose(input);
        write_matrix(matrix_A,n1,n2,stdout);




    }


    MPI_Finalize();

    return 0;
}

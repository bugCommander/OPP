#include <stdio.h>
#include <stdlib.h>
#define a 1e5
#define e 1e-12
#define N 3

#include <math.h>
double phi(double x, double y, double z){
    return x*x+y*y + z*z;
}
double ro(double x,double y,double z){
    return 6 - a*phi(x,y,z);
}


void init_matrix(double *matrix,double *aux_matrix,const double *offset_x,const double *offset_y, const double *offset_z,int size_x,int size_y,int size_z){
    double cur_x,cur_y,cur_z;
    for(int i = 0; i < N;++i){
        cur_z = offset_z[i];
        for(int j = 0; j < N; ++j){
            cur_x = offset_x[j];
            for(int k = 0; k < N; ++k){
                cur_y = offset_y[k];
                int pos = i*N*N + j * N + k;
                if( i == 0 || j ==0 || k == 0 || i == N-1 || j  == N -1 || k == N-1){
                    matrix [ pos] = phi(cur_x,cur_y,cur_z);



                }else {
                    matrix[pos] = 0;
                }
                aux_matrix[pos] = matrix[pos];


            }

        }
    }

}

double max(double x, double y){
    return x>=y?x:y;
}

double
recalk(double *matrix, double *aux_matrix, double *offset_x, double *offset_y, double *offset_z,
       double constant, double delta_x, double delta_y, double delta_z) {

    double max_delta = 0;
    for(int i = 1; i<N-1;++i ){
        for(int j = 1; j < N-1;++j){
            for(int k = 1; k < N-1;++k){
                int pos;
                double x = (aux_matrix[i * N * N + (j - 1) * N + k] + aux_matrix[i * N * N + (j + 1) * N + k]) /
                           (delta_x * delta_x);
                double y = (aux_matrix[i * N * N + j * N + (k - 1)] + aux_matrix[i * N * N + j * N + (k + 1)]) /
                           (delta_y * delta_y);
                double z = (aux_matrix[(i - 1) * N * N + j * N + k] + aux_matrix[(i + 1) * N * N + j * N + k]) /
                           (delta_z * delta_z);
                pos = i * N * N + j * N + k;
                matrix[pos] = constant * (x + y + z - ro(offset_x[i], offset_y[j], offset_z[k]));
                max_delta = max(fabs(matrix[pos] - aux_matrix[pos]), max_delta);


            }
        }
    }
    for(int i = 1;i < N-1; ++i){
        for(int j =1; j < N-1;++j){
            for(int k = 1; k < N-1;++k){
                int pos = i*N*N + j*N + k;
                aux_matrix[pos] = matrix[pos];
            }
        }
    }



    return  max_delta;
    }







double * init_offsets(double min, double delta,int size){
    double *offsets = (double *) malloc(sizeof(double)*size);
    for(int i = 0; i < size;++i){
        offsets[i] = min + delta*i;
    }
    return offsets;

}

double find_big_delta(double *matrix,double *offset_x, double *offset_y, double *offset_z){
    double max_big_delta = 0;
    double aux_big_delta;
    for(int i = 1;i < N-1;++i){
        for(int j =1; j < N-1;++j){
            for(int k = 1; k<N -1;++k){
                aux_big_delta =fabs(matrix[i*N*N+j*N+k] - phi(offset_x[i],offset_y[j],offset_z[k]));
                max_big_delta = max(max_big_delta,aux_big_delta);




            }
        }
    }
    return max_big_delta;


}

int main() {

    const double min_x = -1;
    const double min_y = -1;
    const double min_z = -1;
    const double max_x = 1;
    const double max_y = 1;
    const double max_z =1;
    double delta_x = (max_x - min_x)/(N-1);
    double delta_y = (max_y - min_y)/(N-1);
    double delta_z  = (max_z  - min_z)/(N-1);
    double *offset_x = init_offsets(min_x,delta_x,N);
    double *offset_y = init_offsets(min_y,delta_y,N);
    double *offset_z = init_offsets(min_z,delta_z,N);
    double constant = 1/( 2/(delta_x*delta_x) + 2/(delta_y * delta_y) + 2/(delta_z*delta_z) + a);
    double *matrix = (double*) malloc(sizeof(double)*N*N*N);
    double *aux_matrix = (double*) malloc(sizeof(double)*N*N*N);
    init_matrix(matrix,aux_matrix,offset_x,offset_y,offset_z,N,N,N);
    double max_delta = 0;
    double tmp_delta;

    while(1) {
        tmp_delta = recalk(matrix,aux_matrix,offset_x,offset_y,offset_z,constant,delta_x,delta_y,delta_z);
        if (fabs(tmp_delta - max_delta) < e) {
            break;
        }
        max_delta = tmp_delta;
    }

    double delta = find_big_delta(matrix,offset_x,offset_y,offset_z);
    printf("%lf\n", delta);






/*

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            for(int k = 0; k < N;++k){
                printf("\tmatrix[%d][%d][%d] == %lf \t",j,k,i,matrix[i*N*N + N*j + k]);
            }
            printf("\n");
        }
        printf("\ndafaq\n");

    }
    */



    free(matrix);
    free(aux_matrix);
    free(offset_x);
    free(offset_y);
    free(offset_z);






    return 0;
}

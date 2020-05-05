//
// Created by ksandr on 05.05.2020.
//

#include "inits.h"
void init_proc_list() {
    for (int i = 0; i < size; i++) {
        proc_list[i] = ((i == rank) ? 0 : 1);
    }

}



void init_task_list(int list_size) {
    printf("Creating new task list for rank %d\n", rank);
    list.size = list_size;
    pthread_mutex_lock(&mutex);
    for (int i = 0; i < list.size; i++) {
        list.weights[i] = (rank + 1);
    }
    list.cur_task = 0;
    pthread_mutex_unlock(&mutex);
}
void create_task_list_type() {
    int nitems = 3;
    int blocklengths[] = {1, 1, 1};
    MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(task_list, weights);
    offsets[1] = offsetof(task_list, size);
    offsets[2] = offsetof(task_list, cur_task);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_TASK_LIST);
    MPI_Type_commit(&MPI_TASK_LIST);
}



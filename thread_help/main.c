#include "threads.h"


#include "stddef.h"




int main(int argc, char **argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    create_task_list_type();
    init_task_list(LIST_SIZE);
    proc_list = (int *) malloc(sizeof(int) * size);

    init_proc_list();

    threads_magic();
    free(proc_list);
    MPI_Finalize();
    return 0;
}
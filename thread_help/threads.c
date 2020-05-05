//
// Created by ksandr on 05.05.2020.
//

#include "threads.h"
void do_task(int sig) {
    sleep(sig);

}
int get_extra_task(int target) {
    pthread_mutex_lock(&mutex);
    task_commutator = NEED_TASK;
    pthread_mutex_unlock(&mutex);
    printf(" rank %d  need tasks from rank %d\n", rank, target);
    MPI_Send(&task_commutator, 1, MPI_INT, target, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&task_commutator, 1, MPI_INT, target, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (task_commutator == NO_TASK) {
        printf("rank  %d got no tasks from %d\n", rank, target);
        return NO_TASK;
    }

    printf("rank %d getting tasks   from %d\n", rank, target);
    pthread_mutex_lock(&mutex);
    MPI_Recv(&list, 1, MPI_TASK_LIST, target, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    pthread_mutex_unlock(&mutex);


    return GOT_TASK;
}

int check_recv() {
    for (int i = 0; i < size; i++) {
        if (proc_list[i]) {
            if (get_extra_task(i) == GOT_TASK) {
                return GOT_TASK;
            } else {
                proc_list[i] = 0;
            }
        }
    }

    printf("iteration done for rank %d\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);

    init_task_list(LIST_SIZE);
    init_proc_list();
    return NO_TASK;
}


void *work_thread() {
    printf("Start for iteration %d of %d for rank %d\n", cur_iter + 1, iters, rank);
    while (cur_iter < iters) {
        pthread_mutex_lock(&mutex);
        while (list.cur_task < list.size) {
            printf("Doing task %d of %d with weight %d for rank %d\n",
                   list.cur_task + 1, list.size, list.weights[list.cur_task], rank);
            int task = list.weights[list.cur_task];
            pthread_mutex_unlock(&mutex);
            do_task(task);
            pthread_mutex_lock(&mutex);
            list.cur_task++;
        }
        pthread_mutex_unlock(&mutex);
        printf("Tasks done for rank %d, asking for more\n", rank);
        if (check_recv() == GOT_TASK) {
            continue;
        } else {
            pthread_mutex_lock(&mutex);
            cur_iter++;
            printf("Start for iteration %d of %d for rank %d\n", cur_iter + 1, iters, rank);
            pthread_mutex_unlock(&mutex);
        }
    }
    printf("work done for rank %d\n", rank);
    pthread_mutex_lock(&mutex);
    task_commutator = DEAD;
    pthread_mutex_unlock(&mutex);
    MPI_Send(&task_commutator, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);
    return NULL;
}


void *send_thread() {
    while (cur_iter < iters) {
        MPI_Status status;
        MPI_Recv(&task_commutator, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        printf("im rank %d and i got %d task_commutator\n", rank, task_commutator);
        if (task_commutator == DEAD) {
            printf(" task_commutator DEAD  for rank %d\n", rank);
            break;
        }
        pthread_mutex_lock(&mutex);
        printf("got task_commutator for sending tasks from rank %d, im rank %d\n", status.MPI_SOURCE, rank);
        if (list.cur_task >= list.size - 1) {
            printf("got no tasks to send, im rank %d\n", rank);
            pthread_mutex_unlock(&mutex);
            task_commutator = NO_TASK;
            MPI_Send(&task_commutator, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        } else {
            printf("got tasks to send, im rank %d\n", rank);
            list.size--;
            pthread_mutex_unlock(&mutex);
            task_list new_list;
            new_list.size = 1;
            new_list.cur_task = 0;
            new_list.weights[0] = list.weights[list.size];
            printf("task list from rank %d created\n", rank);
            task_commutator = GOT_TASK;
            printf("sending GOOD task_commutator from rank %d to rank %d \n", rank, status.MPI_SOURCE);
            MPI_Send(&task_commutator, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            printf("sending new tasks from rank %d to rank %d\n", rank, status.MPI_SOURCE);
            MPI_Send(&new_list, 1, MPI_TASK_LIST, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        }


    }
    printf("Data thread rank %d done.\n", rank);
    return NULL;
}


void threads_magic() {
    pthread_mutex_init(&mutex, NULL);
    pthread_attr_t attrs;
    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);
    pthread_t send;
    pthread_t work;
    pthread_create(&send, &attrs, send_thread, NULL);
    pthread_create(&work, &attrs, work_thread, NULL);
    pthread_attr_destroy(&attrs);
    pthread_join(work, NULL);
    pthread_join(send, NULL);


    pthread_mutex_destroy(&mutex);
    MPI_Barrier(MPI_COMM_WORLD);
}
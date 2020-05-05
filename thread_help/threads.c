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
    MPI_Send(&task_commutator, 1, MPI_INT, target, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&task_commutator, 1, MPI_INT, target, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (task_commutator == NO_TASK) {
        return NO_TASK;
    }

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

    MPI_Barrier(MPI_COMM_WORLD);

    init_task_list(LIST_SIZE);
    init_proc_list();
    return NO_TASK;
}


void *work_thread() {
    int aux=0;
    int *dead_weight = malloc(sizeof(int)*size*list.size);
    while (cur_iter < iters) {
        pthread_mutex_lock(&mutex);
        while (list.cur_task < list.size) {
            int task = list.weights[list.cur_task];
            pthread_mutex_unlock(&mutex);
            do_task(task);
            dead_weight[aux++] = task;
        ///   printf("rank %d  ////  weight %d//// number %d\n",rank,task,list.cur_task);
            pthread_mutex_lock(&mutex);
            list.cur_task++;
        }
        pthread_mutex_unlock(&mutex);
        if (check_recv() == GOT_TASK) {
            continue;
        } else {
            pthread_mutex_lock(&mutex);
            printf("rank %d\n",rank);
            for(int i = 0; i < aux;++i){
                printf("solve task /// № %d/// weight %d\n",i, dead_weight[i]);
            }
            cur_iter++;
            aux = 0;
            printf("next iter%d\n",cur_iter);
            pthread_mutex_unlock(&mutex);
        }
    }
    pthread_mutex_lock(&mutex);
    task_commutator = DEAD;
    pthread_mutex_unlock(&mutex);
    MPI_Send(&task_commutator, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);
    free(dead_weight);
    pthread_exit(0);
}


void *send_thread() {
    while (cur_iter < iters) {
        MPI_Status status;
        MPI_Recv(&task_commutator, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        if (task_commutator == DEAD) {
            break;
        }
        pthread_mutex_lock(&mutex);
        if (list.cur_task >= list.size - 1) {
            pthread_mutex_unlock(&mutex);
            task_commutator = NO_TASK;
            MPI_Send(&task_commutator, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        } else {
            list.size--;
            pthread_mutex_unlock(&mutex);
            task_list new_list;
            new_list.size = 1;
            new_list.cur_task = 0;
            new_list.weights[0] = list.weights[list.size-1];
            task_commutator = GOT_TASK;
            MPI_Send(&task_commutator, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(&new_list, 1, MPI_TASK_LIST, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        }


    }
    pthread_exit(0);
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
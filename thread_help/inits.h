//
// Created by ksandr on 05.05.2020.
//

#ifndef THREAD_HELP_INITS_H
#define THREAD_HELP_INITS_H

#include <stdio.h>
#include <mpi/mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#define iters 3
#define REQUEST_TAG 1
#define ANSWER_TAG 2
#define GOT_TASK 3
#define NO_TASK 4
#define NEED_TASK 5
#define DEAD 6
#define MAX_SIZE 50
#define LIST_SIZE 10

typedef struct taskList {
    int weights[MAX_SIZE];
    int size;
    int cur_task;
} task_list;

task_list list;
int size;
int rank;
int task_commutator;
int cur_iter;
MPI_Datatype MPI_TASK_LIST;
int *proc_list;
pthread_mutex_t mutex;


void init_task_list(int list_size);
void init_proc_list();
void create_task_list_type();

#endif //THREAD_HELP_INITS_H

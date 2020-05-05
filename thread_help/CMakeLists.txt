cmake_minimum_required(VERSION 3.15)
project(thread_help C)

set(CMAKE_C_STANDARD 11)

add_executable(thread_help main.c threads.c threads.h inits.c inits.h)

target_link_libraries(thread_help pthread )
target_link_libraries(thread_help mpi)


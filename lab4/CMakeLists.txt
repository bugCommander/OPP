cmake_minimum_required(VERSION 3.15)
project(lab4 C)

set(CMAKE_C_STANDARD 11)

add_executable(non_par non_par.c)
add_executable(dafaq dafaq.c)
target_link_libraries(dafaq mpi)


cmake_minimum_required(VERSION 3.15)
project(lab1 C)

set(CMAKE_C_STANDARD 11)
add_executable(simple simple.c)
add_executable(with_parallel_1 imp1.c)
add_executable(with_parallel_2 imp2.c)

target_link_libraries(simple m)
target_link_libraries(with_parallel_1 m)
target_link_libraries(with_parallel_2 m)



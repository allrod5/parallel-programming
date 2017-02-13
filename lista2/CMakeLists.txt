cmake_minimum_required(VERSION 3.6)
project(lista2)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES ex01.c)

link_libraries(m)

add_executable(lista2 ${SOURCE_FILES})
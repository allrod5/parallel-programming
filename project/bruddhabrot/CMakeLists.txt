cmake_minimum_required(VERSION 3.6)
project(bruddhabrot)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")

set(SOURCE_FILES main.c stb_image_write.h)

link_libraries(m)

add_executable(bruddhabrot ${SOURCE_FILES})
cmake_minimum_required(VERSION 3.6)
project(bruddhabrot)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES main.c)

link_libraries(pthread
               m)

add_executable(bruddhabrot ${SOURCE_FILES})
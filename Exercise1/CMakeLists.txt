cmake_minimum_required(VERSION 3.25)
project(exercise1)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast" )

add_executable(exercise1 exercise1_A.cpp)
cmake_minimum_required(VERSION 3.20)
project(advancedProgramming)

set(CMAKE_CXX_STANDARD 23)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/output )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/output )

add_subdirectory("${PROJECT_SOURCE_DIR}/Exercise1" "${PROJECT_SOURCE_DIR}/Exercise1/Exercise1_output")
add_subdirectory("${PROJECT_SOURCE_DIR}/Exercise2" "${PROJECT_SOURCE_DIR}/Exercise2/Exercise2_output")
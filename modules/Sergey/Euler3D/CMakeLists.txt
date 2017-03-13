cmake_minimum_required(VERSION 3.6)
project(Euler3D)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES src/main.cpp src/Task.cpp src/Task.h src/SparseMatrix.h src/SparseMatrix.cpp)
add_executable(Euler ${SOURCE_FILES})
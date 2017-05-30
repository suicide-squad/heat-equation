cmake_minimum_required(VERSION 2.8)

#find_program(ICC_PATH icc)
#find_program(ICPC_PATH icpc)
#
#set(CMAKE_CXX_COMPILER ${ICPC_PATH})
#set(CMAKE_C_COMPILER ${ICC_PATH})

project(heat-equation)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS}   -std=c11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}    ${OpenMP_C_FLAGS}")
endif()

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

find_library(M_LIB m)

add_subdirectory(modules)

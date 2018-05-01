cmake_minimum_required(VERSION 2.8)

project(heat-equation)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS}   -std=c11")
set(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS}   -D_Float128=__float128")

set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}    ${OpenMP_C_FLAGS}")
endif()

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

# find_library(M_LIB m)

add_subdirectory(modules)

cmake_minimum_required(VERSION 2.8)
project(heat-equation)

if (UNIX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}     -std=c++11")
    set(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS}       -std=c11")
endif(UNIX)


get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE "${DIR_NAME}")

add_subdirectory(modules)
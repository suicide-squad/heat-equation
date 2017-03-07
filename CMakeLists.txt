cmake_minimum_required(VERSION 2.8)
project(heat-equation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE "${DIR_NAME}")

add_subdirectory(modules)
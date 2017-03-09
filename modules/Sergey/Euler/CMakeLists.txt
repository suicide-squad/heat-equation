cmake_minimum_required(VERSION 3.6)
project(Euler)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.cpp src/Task.cpp src/Task.h)
add_executable(Euler ${SOURCE_FILES})
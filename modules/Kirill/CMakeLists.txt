include_directories(include)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE "${DIR_NAME}")

set(MODULE      "${DIR_NAME}")
set(APP         "app_${MODULE}")
set(LIBRARY     "lib_${MODULE}")

set(TYPE        DOUBLE)
set(MKL_RUN     0)
set(AVX2_RUN    1)
set(THREAD_RUN  0)

#   Setting paths
set(INPUT_EULER_SETTING_PATH    "${CMAKE_SOURCE_DIR}/initial/setting3.ini")
set(INPUT_EULER_FUNCTION_PATH   "${CMAKE_SOURCE_DIR}/initial/function3.txt")
set(RESULT_EULER_PATH           "${CMAKE_SOURCE_DIR}/result/Kirill/euler3D_double_avx.txt")
#

find_package(AVX)

if(${AVX2_RUN} AND AVX_FOUND)
    set(AVX2_RUN, 1)
    message(STATUS "AVX2 RUN")
else()
    set(AVX2_RUN, 0)
endif()

find_package(MKL)
if(${MKL_RUN} AND ${MKL_FOUND})
    set(MKL_RUN, 1)
    message(STATUS "MKL RUN")
else()
    set(MKL_RUN, 0)
endif()

configure_file(
        "${CMAKE_CURRENT_LIST_DIR}/include/ts.h.in"
        "${CMAKE_CURRENT_LIST_DIR}/include/ts.h"
)

add_subdirectory(src)
add_subdirectory(app)
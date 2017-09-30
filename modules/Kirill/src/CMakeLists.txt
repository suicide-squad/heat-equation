set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")


if(${MKL_RUN})
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
endif()
if(${AVX2_RUN})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
endif()
add_library(${target} STATIC ${hrcs} ${srcs})
target_link_libraries(${target} ${MPI_C_LIBRARIES})
target_link_libraries(${target} ${M_LIB})
target_link_libraries(${target} ${MKL_LIBRARIES})
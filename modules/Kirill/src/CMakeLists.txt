set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

find_package(MKL)

if (MKL_FOUND AND ${MKL_RUN})
    add_definitions(-DMKL_RUN)
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
endif ()

#add_definitions(-DDOUBLE_TYPE)

add_library(${target} STATIC ${hrcs} ${srcs})
target_link_libraries(${target} ${MPI_C_LIBRARIES})
target_link_libraries(${target} ${M_LIB})
target_link_libraries(${target} ${MKL_LIBRARIES})
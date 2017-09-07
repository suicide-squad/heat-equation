set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

add_library(${target} STATIC ${hrcs} ${srcs})
target_link_libraries(${target} ${MPI_C_LIBRARIES})
target_link_libraries(${target} ${M_LIB})
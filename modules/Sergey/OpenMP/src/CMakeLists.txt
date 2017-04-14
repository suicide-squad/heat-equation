set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(${target} STATIC ${hrcs} ${srcs})

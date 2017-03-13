set(target ${APP})

file(GLOB  src "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

add_executable(${target} ${src})
target_link_libraries(${target} ${LIBRARY})

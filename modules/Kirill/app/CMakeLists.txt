set(target ${APP})

file( GLOB SOURCE_DOUBLE_EULER "main.c" )
add_executable(${target} ${SOURCE_DOUBLE_EULER})
target_link_libraries(${target} ${LIBRARY})

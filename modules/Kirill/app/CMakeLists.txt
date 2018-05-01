# ЯВНЫЙ МЕТОД
set(target "${APP}_euler")

file(GLOB  src_euler euler.c)
include_directories(${OpenCL_INCLUDE_DIRS})
add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})

## РУНГЕ-КУТТ
#set(target "${APP}_runge")
#
#file(GLOB  src_runge runge.c)
#
#add_executable(${target} ${src_runge})
#target_link_libraries(${target} ${LIBRARY})
#target_link_libraries(${target} ${MPI_C_LIBRARIES})
#
### НЕЯВНЫЙ МЕТОД
#set(target "${APP}_implicit")
#
#file(GLOB  src_implicit implicit.c)
#
#add_executable(${target} ${src_implicit})
#target_link_libraries(${target} ${LIBRARY})

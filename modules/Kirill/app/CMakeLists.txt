# ЯВНЫЙ МЕТОД
set(target "${APP}_euler")

file(GLOB  src_euler euler.c)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_C_LIBRARIES})

# НАИВНЫЙ ЯВНЫЙ МЕТОД
set(target "${APP}_euler_naive")

file(GLOB  src_euler_naive euler_naive.c)

add_executable(${target} ${src_euler_naive})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_C_LIBRARIES})

## РУНГЕ-КУТТ
set(target "${APP}_runge")

file(GLOB  src_runge runge.c)

add_executable(${target} ${src_runge})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_C_LIBRARIES})

## НАИВНЫЙ РУНГЕ-КУТТ
set(target "${APP}_runge_naive")

file(GLOB  src_runge_naive runge_naive.c)

add_executable(${target} ${src_runge_naive})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_C_LIBRARIES})

## НЕЯВНЫЙ МЕТОД
set(target "${APP}_implicit")

file(GLOB  src_implicit implicit.c)

add_executable(${target} ${src_implicit})
target_link_libraries(${target} ${MPI_C_LIBRARIES})
target_link_libraries(${target} ${LIBRARY})

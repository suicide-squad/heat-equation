# MPI TEST
set(target "${APP}_Test")

file(GLOB src_mpi_test mpi_test.cpp)

add_executable(${target} ${src_mpi_test})
target_link_libraries(${target} ${MPI_LIBRARIES})

# ЯВНЫЙ МЕТОД
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})



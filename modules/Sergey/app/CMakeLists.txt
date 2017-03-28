# ЯВНЫЙ МЕТОД
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})

# Рунхе Хутта
set(target "${APP}_Runge")

file(GLOB src_runge Runge.cpp)

add_executable(${target} ${src_runge})
target_link_libraries(${target} ${LIBRARY})
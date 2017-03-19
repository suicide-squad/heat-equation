# ЯВНЫЙ МЕТОД
set(target "${APP}_euler")

file(GLOB  src_euler euler.c)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})


# РУНГЕ-КУТТ
set(target "${APP}_runge")

file(GLOB  src_runge runge.c)

add_executable(${target} ${src_runge})
target_link_libraries(${target} ${LIBRARY})

## НЕЯВНЫЙ МЕТОД
#set(target "${APP}_imp")
#
#file(GLOB  src_imp imp.c)
#
#add_executable(${target} ${src_imp})
#target_link_libraries(${target} ${LIBRARY})

set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

add_library(${target} STATIC ${hrcs} ${srcs})

target_link_libraries(${target} ${MPI_C_LIBRARIES})
target_link_libraries(${target} ${M_LIB})

if(${MKL_RUN})
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
    target_link_libraries(${target} ${MKL_LIBRARIES})
elseif(${AVX2_RUN})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
elseif(${CPUGPU_RUN})
	include_directories(${OpenCL_INCLUDE_DIRS})
    target_link_libraries(${target} ${OPENCL_LIBRARIES})
elseif(${FPGA_RUN})
	include_directories(${OpenCL_INCLUDE_DIRS})
    target_link_libraries(${target} ${OPENCL_LIBRARIES})
   	target_link_libraries(${target} ${AlteraCL_LIBRARY})
endif()

file(GLOB src_cl "*.cl")
file(COPY ${src_cl} DESTINATION ".")
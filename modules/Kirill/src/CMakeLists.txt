set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

if (${CPU_CL_RUN} OR ${GPU_CL_RUN} OR ${FPGA_RUN})
    list(APPEND srcs "${CMAKE_CURRENT_SOURCE_DIR}/utils/opencl.c")
endif()

message(STATUS ${srcs})

add_library(${target} STATIC ${hrcs} ${srcs})

target_link_libraries(${target} ${M_LIB})

# Added MPI Libraries
if (NOT ${LIB} MATCHES "opencl")
    message(STATUS "MPI LIB")
else()
    target_link_libraries(${target} ${MPI_C_LIBRARIES})
endif()

if(${MKL_RUN})
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
    target_link_libraries(${target} ${MKL_LIBRARIES})
elseif(${AVX2_RUN})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native ")
elseif(${CPU_CL_RUN} OR ${GPU_CL_RUN})
    message(STATUS "OPENCL_INCLUDE_DIRS - ${OPENCL_INCLUDE_DIRS}")
    include_directories(${OPENCL_INCLUDE_DIRS})
    target_link_libraries(${target} ${OPENCL_LIBRARIES})
elseif(${FPGA_RUN})
	include_directories(${OpenCL_INCLUDE_DIRS})
    target_link_libraries(${target} ${OPENCL_LIBRARIES})
   	target_link_libraries(${target} ${AlteraCL_LIBRARY})
endif()

# Copy file kernel
file(GLOB src_cl "*.cl")
file(COPY ${src_cl} DESTINATION ".")

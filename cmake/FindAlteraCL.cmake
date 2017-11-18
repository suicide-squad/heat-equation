FIND_PATH(OpenCL_INCLUDE_DIRS
        NAMES
        CL/cl.h CL/opencl.h
        PATHS
        $ENV{ALTERAOCLSDKROOT}/host/include
        ${ALTERAOCLSDKROOT}/host/include
        $ENV{${ALTERAOCLSDKROOT}}/host/include)

FIND_LIBRARY(OpenCL_LIBRARIES
        NAMES
        OpenCL
        PATHS
        $ENV{ALTERAOCLSDKROOT}/linux64/lib
        ${ALTERAOCLSDKROOT}/linux64/lib
        $ENV{${ALTERAOCLSDKROOT}}/linux64/lib)

FIND_LIBRARY(AlteraCL_LIBRARY
        NAMES
        libalteracl.so
        PATHS
        $ENV{ALTERAOCLSDKROOT}/linux64/lib
        ${ALTERAOCLSDKROOT}/linux64/lib
        $ENV{${ALTERAOCLSDKROOT}}/linux64/lib)

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
        $ENV{ALTERAOCLSDKROOT}/host/arm32/lib
        ${ALTERAOCLSDKROOT}/host/arm32/lib
        $ENV{${ALTERAOCLSDKROOT}}/host/arm32/lib)

FIND_LIBRARY(AlteraCL_LIBRARY
        NAMES
        libalteracl.so
        PATHS
        $ENV{ALTERAOCLSDKROOT}/board/de10_nano/arm32/lib
        ${ALTERAOCLSDKROOT}/board/de10_nano/arm32/lib
        $ENV{${ALTERAOCLSDKROOT}}/board/de10_nano/arm32/lib)

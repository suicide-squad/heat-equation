
FIND_LIBRARY(AlteraCL_LIBRARY
        NAMES
        libalteracl.so
        PATHS
        $ENV{ALTERAOCLSDKROOT}/linux64/lib
        ${ALTERAOCLSDKROOT}/linux64/lib
        $ENV{${ALTERAOCLSDKROOT}}/linux64/lib)

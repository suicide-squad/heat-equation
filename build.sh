#!/bin/bash

Clear() {
  if [ -d _build ]; then
    rm -rf _build
  fi
}

while [[ $# -ge 1 ]]
do
    key=$1
    case $key in
        -h|--help)
            echo "Flags"
            echo "-a|--arc"
            echo "-c|--compiler"
            echo "-l|--lib"
            echo "-f|--filename"
            echo "FPGA compile $./build.sh -a fpga -l opencl -f result.txt"
            exit
            ;;
        -c|--compiler)
            COMPILER=$2
            shift
            ;;
        -a|--arc)
            ARCH=$2
            shift
            ;;
        -l|--lib)
            LIB=$2
            shift
            ;;
        -f|--filename)
            OUT_FILE=$2
            shift
            ;;
        -m|--mode)
            BUILD_MODE=$2
            shift
            ;;
        -r|--rm)
            Clear
            shift
            ;;
        -ca|--compile-aocx)
            COMPILE_AOCX_FLAG=$2
            shift
            ;;
        *)
            ;;
    esac

    shift
done

case $COMPILER in
    intel)
        source /opt/intel/bin/compilervars.sh intel64
        C_COMPILER="icc"
        CXX_COMPILER="icpc"
        ;;
    *)
    # TO-DO Get back a configuration for emulator. Also need revert in BuildAltera
        if [ "$ARCH" = "fpga" ]; then
            C_COMPILER="arm-linux-gnueabihf-gcc"
            CXX_COMPILER="arm-linux-gnueabihf-g++"
        else
            C_COMPILER="gcc"
            CXX_COMPILER="g++"
        fi
        ;;
esac

if [ "$ARCH" == "gpu" ] || [ "$ARCH" == "fpga" ]
then
    LIB=opencl
elif [ "$ARCH" == "cpu" ] && [ "$LIB" == "mkl" ] && [ "$C_COMPILER" == "icc" ]
then
    LIB=mkl
elif [ "$ARCH" == "cpu" ] && [ "$LIB" == "omp" ]
then
    LIB=omp
elif [ "$ARCH" == "cpu" ] && [ "$LIB" == "avx" ]
then
    LIB=avx
fi

if [ "$LIB" == "" ]; then
    LIB=""
fi

if [ "$BUILD_MODE" == "debug" ]; then
    BUILD_MODE="Debug"
else
    BUILD_MODE="Release"
fi

if [ "$OUT_FILE" == "" ]
then
    OUT_FILE=euler3D.txt
fi


Build() {
    if [ "$ARCH" = "fpga" ]; then
        echo "Run with fpga"
        source env_altera.sh
    fi

    root_dir="$( pwd )"
    if [ ! -d _build ]; then
        mkdir _build
    fi
    cd _build
    if [ -f "CMakeCache.txt" ]; then
        rm "CMakeCache.txt"
    fi

    echo "C_COMPILER = $C_COMPILER"
    echo "CXX_COMPILER = $CXX_COMPILER"
    cmake -DCMAKE_BUILD_TYPE=$BUILD_MODE -DCMAKE_C_COMPILER=$C_COMPILER -DCMAKE_CXX_COMPILER=$CXX_COMPILER -DARCH=$ARCH -DLIB=$LIB -DOUT_FILE=$OUT_FILE  $root_dir
    make -j2
    cd $root_dir
}

BuildAltera() {
    cd _build/modules/Kirill/src
    aoc kernel.cl -o ./kernel.aocx -board=de10_nano_sharedonly --profile
#    CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1 ./Hello
}

Build
if [ "$ARCH" = "fpga" ] && [ "$COMPILE_AOCX_FLAG" = 1 ]; then
    BuildAltera
fi

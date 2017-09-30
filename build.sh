#!/bin/bash

COMPILER=$1
BUILD=$2

Build() {
    C_COMPILER="gcc"
    CXX_COMPILER="g++"
    TYPE_BUILD=Release
    root_dir="$( pwd )"
    if [ ! -d _build ]; then
        mkdir _build
    fi
    cd _build
    if [ -f "CMakeCache.txt" ]; then
        rm "CMakeCache.txt"
    fi
    echo $root_dir
    if [ "$COMPILER" = "intel" ]; then
      source /opt/intel/bin/iccvars.sh intel64
      C_COMPILER="icc"
      CXX_COMPILER="icpc"
    fi
    if [ "$BUILD" = "debug" ]; then
        TYPE_BUILD=Debug
    fi
    echo "COMPILER = $C_COMPILER"
    echo "TYPE_BUILD = $TYPE_BUILD"
    cmake -DCMAKE_BUILD_TYPE=$TYPE_BUILD -DCMAKE_C_COMPILER=$C_COMPILER -DCMAKE_CXX_COMPILER=$CXX_COMPILER $root_dir
    make -j2
    cd $root_dir
}

Clear() {
  if [ -d _build ]; then
    rm -rf _build
  fi
}

if [ "$3" = "clear" ]; then
    echo "run clear directory"
    Clear
#    echo "finish clear directory"
fi

Build

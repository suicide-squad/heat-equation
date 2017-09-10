#!/bin/bash

Build() {
    C_COMPILER="gcc"
    CXX_COMPILER="g++"

    root_dir="$( pwd )"
    if [ ! -d build ]; then
        mkdir build
    fi
    cd build
    echo $root_dir
    if [ "$1" = "intel" ]; then
      source /opt/intel/bin/iccvars.sh intel64
      C_COMPILER="icc"
      CXX_COMPILER="icpc"
    fi
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=$C_COMPILER -DCMAKE_CXX_COMPILER=$CXX_COMPILER $root_dir
    make -j2
    cd $root_dir
}

Clear() {
  if [ -d build ]; then
    rm -rf ./build
  fi
}

if [ "$2"="clear" ]; then
    Clear
fi

Build

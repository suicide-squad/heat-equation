#!/bin/bash

COMPILER=$1

Build() {
    C_COMPILER="gcc"
    CXX_COMPILER="g++"

    root_dir="$( pwd )"
    if [ ! -d _build ]; then
        mkdir _build
    fi
    cd _build
    echo $root_dir
    if [ "$COMPILER" = "intel" ]; then
      source /opt/intel/bin/iccvars.sh intel64
      C_COMPILER="icc"
      CXX_COMPILER="icpc"
    fi
    echo "COMPILER = $C_COMPILER"
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=$C_COMPILER -DCMAKE_CXX_COMPILER=$CXX_COMPILER $root_dir
    make -j2
    cd $root_dir
}

Clear() {
  if [ -d _build ]; then
    rm -rf _build
  fi
}

if [ "$2" = "clear" ]; then
    echo "run clear directory"
    Clear
#    echo "finish clear directory"
fi

Build

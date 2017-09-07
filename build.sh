#!/bin/sh

Build () {
    root_dir="$( pwd )"
    mkdir build
    cd build
    echo $root_dir
    cmake -DCMAKE_BUILD_TYPE=RELEASE $root_dir
    make -j2
}

Clear() {
  if [ -d build ]
    then
    rm -rf ./build
  fi
}

Clear
Build

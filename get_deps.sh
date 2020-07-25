#!/bin/bash
set -e -x

if [ ! -d exaCore ]; then
  git clone https://github.com/exaAlgo/exaCore.git
fi

if [ ! -d exaSort ]; then
  git clone https://github.com/exaAlgo/exaSort.git -b refactor
fi

cd exaCore
make MPI="${MPI}" GS_DIR="${GS_DIR}" OCCA=0 PREFIX="${PREFIX}" install

export EXA_DEBUG=0
cd ../exaSort
make MPI="${MPI}" GS_DIR="${GS_DIR}" EXA_DIR="${PREFIX}" PREFIX="${PREFIX}" install
make MPI="${MPI}" GS_DIR="${GS_DIR}" EXA_DIR="${PREFIX}" PREFIX="${PREFIX}" tests

set +x

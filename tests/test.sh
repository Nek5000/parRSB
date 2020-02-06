#!/bin/bash

export CASE=bp5
export VERBOSE=0

export NEK_SOURCE_ROOT=`pwd`/Nek5000
export PPLIST="PARRSB DPROCMAP"

export LD_LIBRARY_PATH="${BUILDDIR}:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${EXADIR}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${EXASORTDIR}/lib:${LD_LIBRARY_PATH}"

function print_test_err(){
  echo "Test: $1, not ok."
}

function print_test_success(){
  echo "Test: $1, ok."
}

function exit_err(){
  echo $1
  exit 1
}

function clean_nek(){
  rm -rf *.co2 *.f build.log makefile nek5000\
    obj/ SESSION.NAME 2>/dev/null *.out *.err 2>/dev/null
}

function clean_con(){
  rm *.co2 parCon *.out *.err *.run 2>/dev/null
}

function gen_box(){
  ${NEK_SOURCE_ROOT}/bin/genbox >log << EOF
genbox.in
EOF
  mv box.re2 ${CASE}.re2
}

function setup_test(){
  cp ./bbbb/* ${NEK_SOURCE_ROOT}/bin/makenek $1/;

  cd $1;
  gen_box;
  cd - 2>/dev/null 1>&2
}

function setup_nek5000(){
  # TODO: Copy a specific nek release
  git clone https://github.com/Nek5000/Nek5000.git

  cd ${NEK_SOURCE_ROOT}/tools
  ./maketools genbox gencon
  cd - 2>/dev/null 1>&2
}

function run_nek5000(){
  local CASE=$1
  local NP=$2

  echo ${CASE}   >  SESSION.NAME
  echo `pwd`'/' >>  SESSION.NAME
  mpirun -np ${NP} ./nek5000 >nek.out 2>nek.err
}

function run_test(){
  local TEST=$1
  local NP=$2
  local ORIGDIR=`pwd`

  cd ${TEST}

  if [ ${VERBOSE} -eq 1 ]; then
    echo "Building ${TEST} ..."
  fi
  ./makenek ${CASE} >build.out 2>build.err
  if [ -s build.err ]; then
    cd ${ORIGDIR}
    exit_err "Error building ${TEST}."
  fi

  if [ ${VERBOSE} -eq 1 ]; then
    echo "Copying parCon from ${BUILDDIR}/examples/parCon"
  fi
  cp ${BUILDDIR}/examples/parCon . 2>/dev/null
  if [ ! -f ./parCon ]; then
    cd ${ORIGDIR}
    exit_err "Error copying parCon."
  fi

  if [ ${VERBOSE} -eq 1 ]; then
    echo "Running parCon on ${CASE}.re2 ..."
  fi
  if [ ${NP} -gt 8 ]; then
    NP=8
  fi
  mpirun -np ${NP} ./parCon ${CASE}.re2 \
    >parcon.out 2>parcon.err
  if [ -s parcon.err ]; then
    cd ${ORIGDIR}
    exit_err "Error running parCon."
  fi

  if [ ${VERBOSE} -eq 1 ]; then
    echo "Running ${TEST} ..."
  fi
  run_nek5000 ${CASE} ${NP}

  success=`grep 'run successful: dying ...' nek.out | wc -l`
  if [ ${success} -eq 0 ]; then
    print_test_err "${TEST}"
  else
    print_test_success "${TEST}"
  fi

  cd ${ORIGDIR}
}

function run_test_suite(){
  for i in ./[bp]0*; do
    setup_test $i

    NP=${i:(-3)};
    run_test ${i} ${NP}
  done
}

while (( "$#" )); do
  case "$1" in
    -v|--verbose)
      export VERBOSE=1
      shift
      ;;
    --get-deps)
      setup_nek5000
      shift
      ;;
    -r|--run)
      run_test_suite
      shift
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

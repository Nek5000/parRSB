#!/bin/bash

export CASE=nek-case
export VERBOSE=0

export NEK_SOURCE_ROOT=`pwd`/Nek5000
export PPLIST="PARRSB DPROCMAP"

export LD_LIBRARY_PATH="${BUILDDIR}:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${EXADIR}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${EXASORTDIR}/lib:${LD_LIBRARY_PATH}"

dir_resolve()
{
  cd "$1" 2>/dev/null || return $?  # cd to desired directory;
  #if fail, quell any error messages but return exit status
  echo "`pwd -P`" # output full, link-resolved path
}

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

function setup_nek5000(){
  # TODO: Copy a specific nek release
  git clone https://github.com/Nek5000/Nek5000.git

  cd ${NEK_SOURCE_ROOT}/tools
  ./maketools genbox gencon
  cd - 2>/dev/null 1>&2
}

function clean_con(){
  rm *.co2 *.re2 2>/dev/null
}

function gen_box(){
  ${NEK_SOURCE_ROOT}/bin/genbox >log << EOF
genbox.in
EOF
  mv box.re2 ${CASE}.re2
}

function setup_test(){
  local TEST=$1

  clean_con

  cp ${BUILDDIR}/examples/parCon . 2>/dev/null
  if [ ! -f ./parCon ]; then
    exit_err "Error copying parCon."
  fi

  cp ${NEK_SOURCE_ROOT}/bin/makenek ${TEST}/genbox.in .
  gen_box

  ./makenek ${CASE} >build.out 2>build.err
  if [ -s build.err ]; then
    exit_err "Error building ${TEST}."
  fi
}

function run_nek5000(){
  local NP=$1
  local TEST=$2

  echo ${CASE}   >  SESSION.NAME
  echo `pwd`'/' >>  SESSION.NAME
  mpirun -np ${NP} ./nek5000 >nek.${TEST}.out 2>nek.${TEST}.err
}

function run_test(){
  local TEST=$1
  local NP=1

  if [ $# -gt 1 ]; then NP=$2;
  else NP=${i:(-3)}; fi

  setup_test $i

  if [ ${NP} -gt 8 ]; then NP=7; fi
  mpirun -np ${NP} ./parCon ${CASE}.re2 >parcon.out 2>parcon.err
  if [ -s parcon.err ]; then
    exit_err "Error running parCon."
  fi

  TEST="`dir_resolve \"${TEST}\"`"
  TEST=`basename ${TEST}`
  run_nek5000 ${NP} ${TEST}

  success=`grep 'run successful: dying ' nek.${TEST}.out | wc -l`
  if [ ${success} -eq 0 ]; then
    print_test_err "${TEST}"
  else
    print_test_success "${TEST}"
  fi
}

function run_test_suite(){
  cd ${CASE}
  # Run 2d tests first
  cp SIZE.2d SIZE
  for i in ../[wp]2d*; do
    run_test ${i}
  done
  # Run 3d tests
  cp SIZE.3d SIZE
  for i in ../[wp]3d*; do
    run_test ${i}
  done
  cd - 2>/dev/null 1>&2
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

#!/bin/bash

export VERBOSE=0

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

  export NEK_SOURCE_ROOT=`pwd`/Nek5000

  cd ${NEK_SOURCE_ROOT}/tools
  ./maketools genbox gencon
  cd - 2>/dev/null 1>&2
}

function gen_box(){
  local test=$1

  ${NEK_SOURCE_ROOT}/bin/genbox >log << EOF
genbox.in
EOF
}

function gen_con(){
  local rea=$1

  ${NEK_SOURCE_ROOT}/bin/gencon >log << EOF
${rea}
0.01
EOF
}

function run_test(){
  local test=$1
  local np=1

  # Set number of processors for the test
  if [ $# -gt 1 ]; then np=$2;
  else np=${i:(-3)}; fi
  if [ ${np} -gt 8 ]; then np=7; fi

  cd ${test}

  # clean .co2 or .re2 files if present
  rm *.co2 *.re2 2>/dev/null

  # generate the geometry
  gen_box ${test}
  mv box.re2 ${test}.re2

  #generate the .co2 file by gencon (Remove this at some point)
  gen_con ${test}
  mv ${test}.co2 gencon.co2

  # copy parCOn
  cp ${BUILDDIR}/examples/parCon . 2>/dev/null
  if [ ! -f ./parCon ]; then
    exit_err "Error copying parCon."
  fi

  # run parCon to generate the .co2 file
  mpirun -np ${np} ./parCon ./${test}.re2 >out 2>err
  mv ${test}.co2 parCon.co2

  hash1=(`md5sum gencon.co2`)
  hash2=(`md5sum parCon.co2`)

  if [ "${hash1}" = "${hash2}" ]; then
    print_test_success "${test}"
  else
    print_test_err "${test}"
  fi

  cd - 2>&1 > /dev/null
}

function run_test_suite(){
  setup_nek5000

  for i in [wp]*; do
    run_test ${i}
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

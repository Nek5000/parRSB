#!/bin/bash

export VERBOSE=0

function print_test_failure(){
  local RED='\033[0;31m'
  local NC='\033[0m'
  echo -e "${RED} $1 ${NC}"
}

function print_test_success(){
  local GREEN='\033[0;92m'
  local NC='\033[0m'
  echo -e "${GREEN} $1 ${NC}"
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

function gen_mesh(){
  local test=$1

  # clean .co2 or .re2 files if present
  rm *.co2 *.re2 2>/dev/null

  # generate the geometry
  ${NEK_SOURCE_ROOT}/bin/genbox >log << EOF
genbox.in
EOF
  mv box.re2 ${test}.re2

  #generate the .co2 file by gencon (Remove this at some point)
  ${NEK_SOURCE_ROOT}/bin/gencon >log << EOF
${test}
0.01
EOF

}

function run_test_suite(){
  setup_nek5000

  tests=(`ls t[0-9][0-9][0-9]*`)
  meshes=(`ls -d w[23]d_*`)

  echo "${tests[@]}"
  echo "${meshes[@]}"

  for m in "${meshes[@]}"; do
    cd ${m}
    gen_mesh ${m}
    cd - 2>&1 >/dev/null
  done

  for t in "${tests[@]}"; do
    for m in "${meshes[@]}"; do
      cd ${m}
      mpirun -np 1 ../${t} "`pwd`/${m}.co2" >out.log 2>err.log
      wait $!
      if [ ! -s err.log ]; then
        print_test_success "Test: ${t} ${b} ok."
      else
        print_test_failure "Test: ${t} ${b} not ok."
      fi
      cd - 2>&1 >/dev/null
    done
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

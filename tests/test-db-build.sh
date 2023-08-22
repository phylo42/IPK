#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "Parameters: $#"

if [ $# -eq 2 ]; then
    echo "Parameters: $@"
    ROOT_DIR=`realpath $1`
    RAXML_NG=`realpath $2`
else
    ROOT_DIR=`realpath "${SCRIPT_DIR}"/../..`
    RAXML_NG=`which raxml-ng`
fi



BIN_DIR="${ROOT_DIR}"/bin
IPK_BIN="${BIN_DIR}"/ipk/ipk-dna
IPK_SCRIPT="${ROOT_DIR}"/ipk.py
IPK_DIFF_BIN="${BIN_DIR}"/tools/ipkdiff-dna
WORKING_DIR="${ROOT_DIR}"/output

echo "Pwd: `pwd`"
echo "Root dir: ${ROOT_DIR}"
echo "Bin dir: ${BIN_DIR}"
echo



# Neotrop test
NEOTROP_REFERENCE="${SCRIPT_DIR}"/data/neotrop/reference.fasta
NEOTROP_TREE="${SCRIPT_DIR}"/data/neotrop/tree.rooted.newick
NEOTROP_DATABASE_REFERENCE="${SCRIPT_DIR}"/data/neotrop/DB_k7_o2.0.ipk

# D140 test (proteins)
D140_REFERENCE="${SCRIPT_DIR}"/data/D140/reference.fasta
D140_TREE="${SCRIPT_DIR}"/data/D140/tree.newick
D140_DATABASE_REFERENCE="${SCRIPT_DIR}"/data/D140/DB_k4_o10.ipk


function run_test() {
    local REFERENCE=$1
    local TREE=$2
    local DATABASE_REFERENCE=$3
    local DATABASE_BUILD="${WORKING_DIR}"/$(basename "${DATABASE_REFERENCE}")

    mkdir -p "${WORKING_DIR}"
    rm -f "${DATABASE_BUILD}"

    local command=python3 "${IPK_SCRIPT}" build -r "${REFERENCE}" -t "${TREE}" -m GTR -k 7 --omega 2.0 -u 1.0 -b "${RAXML_NG}" -w "${WORKING_DIR}"

    echo "Binary files: OK. Running IPK as: ${command}"
    eval "${command}"

    if [ ! -f "${DATABASE_BUILD}" ]; then
        echo "Error: could not find ${DATABASE_BUILD}. Something went wrong"
        exit 4
    fi

    $IPK_DIFF_BIN 0 "${DATABASE_REFERENCE}" "${DATABASE_BUILD}"

    if [ $? -ne 0 ]; then
        echo "Error: databases are different. See the ipkdiff log"
        exit 5
    fi
}

# Check binary files and tools
if [ ! -f "${IPK_BIN}" ]; then
    echo "Error: could not find binary files of IPK: ${IPK_BIN}. Please make sure to compile the project"
    exit 1
elif [ ! -f "${IPK_DIFF_BIN}" ]; then
    echo "Error: could not find tools: ${IPK_DIFF_BIN}. Please make sure to compile it separately, i.e. do 'make diff-dna' or 'cmake --build DIR --target diff-dna"
    exit 2
elif [ ! "${RAXML_NG}" ]; then
    echo "Error: could not find raxml-ng."
    exit 3
else
    # Run tests
    run_test "${NEOTROP_REFERENCE}" "${NEOTROP_TREE}" "${NEOTROP_DATABASE_REFERENCE}"
    run_test "${ANOTHER_REFERENCE}" "${ANOTHER_TREE}" "${ANOTHER_DATABASE_REFERENCE}"
fi


if [ ! -f "${IPK_BIN}" ]
then
    echo "Error: could not find binary files of IPK: ${IPK_BIN}. Please make sure to compile the project"
    exit 1
elif [ ! -f "${IPK_DIFF_BIN}" ]
then
    echo "Error: could not find tools: ${IPK_DIFF_BIN}. Please make sure to compile it separately, i.e. do 'make diff-dna' or 'cmake --build DIR --target diff-dna"
    exit 2
elif [ ! "${RAXML_NG}" ]
then
    echo "Error: could not find raxml-ng."
    exit 3
else

    # Run tests
    run_test "${NEOTROP_REFERENCE}" "${NEOTROP_TREE}" "${NEOTROP_DATABASE_REFERENCE}"
    run_test "${D140_REFERENCE}" "${D140_TREE}" "${D140_DATABASE_REFERENCE}"
fi
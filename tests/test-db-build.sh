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

exit
REFERENCE="${SCRIPT_DIR}"/data/neotrop/reference.fasta
TREE="${SCRIPT_DIR}"/data/neotrop/tree.rooted.newick
DATABASE_REFERENCE="${SCRIPT_DIR}"/data/neotrop/DB_k7_o2.0.rps
DATABASE_BUILD="${WORKING_DIR}"/DB_k7_o2.0.rps


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
    mkdir -p "${WORKING_DIR}"
    rm -f "${DATABASE_BUILD}"


    command=python3 "${IPK_SCRIPT}" build -r "${REFERENCE}" -t "${TREE}" -m GTR -k 7 --omega 2.0 -u 1.0 -b "${RAXML_NG}" -w "${WORKING_DIR}"

    echo "Binary files: OK. Running IPK as: ${command}"
    eval "${command}"

    if [ ! -f "${DATABASE_BUILD}" ]
    then
        echo "Error: could not find ${DATABASE_BUILD}. Something went wrong"
        exit 4
    fi

    $IPK_DIFF_BIN 0 "${DATABASE_REFERENCE}" "${DATABASE_BUILD}"

    if [ $? -ne 0 ]; then
        echo "Error: databases are different. See the ipkdiff log"
        exit 5
    fi
fi
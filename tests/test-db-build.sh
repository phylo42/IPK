#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "Parameters: $#"

if [ $# -eq 2 ]; then
    echo "Parameters: $@"
    ROOT_DIR=`realpath $1`
    WORKING_DIR=`realpath $2`
    RAXML_NG=$ROOT_DIR/raxml-ng
    BIN_DIR="${ROOT_DIR}"
else
    ROOT_DIR=`realpath "${SCRIPT_DIR}"/..`
    RAXML_NG=`which raxml-ng`
    WORKING_DIR="${ROOT_DIR}"/output
fi


BIN_DIR="${ROOT_DIR}"
IPK_BIN="${BIN_DIR}"/ipk-dna
IPK_SCRIPT="${ROOT_DIR}"/ipk.py
IPK_DIFF_BIN="${BIN_DIR}"/ipkdiff-dna
IPK_DIFF_AA_BIN="${BIN_DIR}"/ipkdiff-aa

echo "Pwd: `pwd`"
echo "Root dir: ${ROOT_DIR}"
echo "Bin dir: ${BIN_DIR}"
echo


if [ ! -f "${IPK_BIN}" ]
then
    echo "Error: could not find binary files of IPK: ${IPK_BIN}. Please compile and install IPK"
    exit 1
elif [ ! -f "${IPK_SCRIPT}" ]
then
    echo "Error: could not find ${IPK_SCRIPT}. Something went wrong"
    exit 2
elif [ ! -f "${IPK_DIFF_BIN}" ]
then
    echo "Error: could not find tools: ${IPK_DIFF_BIN}. Please make sure to compile it separately, i.e. do 'make diff-dna' or 'cmake --build DIR --target diff-dna"
    exit 3
elif [ ! "${RAXML_NG}" ]
then
    echo "Error: could not find raxml-ng."
    exit 4
else
    mkdir -p "${WORKING_DIR}"

    # Neotrop test
    NEOTROP_REFERENCE="${SCRIPT_DIR}"/data/neotrop/reference.fasta
    NEOTROP_TREE="${SCRIPT_DIR}"/data/neotrop/tree.rooted.newick
    NEOTROP_DATABASE_REFERENCE="${SCRIPT_DIR}"/data/neotrop/DB_k7_o2.0.ipk
    NEOTROP_DATABASE_BUILD="${WORKING_DIR}"/DB_k7_o2.0.ipk

    rm -f "${NEOTROP_DATABASE_BUILD}"
    command="${IPK_SCRIPT}" build -r "${NEOTROP_REFERENCE}" -t "${NEOTROP_TREE}" -m GTR -k 7 --omega 2.0 -b "${RAXML_NG}" -w "${WORKING_DIR}" -o "${NEOTROP_DATABASE_BUILD}"

    echo "Binary files: OK. Running IPK as: ${command}"
    eval "${command}"

    if [ ! -f "${NEOTROP_DATABASE_BUILD}" ]
    then
        echo "Error: could not find ${NEOTROP_DATABASE_BUILD}. Something went wrong"
        exit 5
    fi

    $IPK_DIFF_BIN 0 "${NEOTROP_DATABASE_REFERENCE}" "${NEOTROP_DATABASE_BUILD}"

    if [ $? -ne 0 ]; then
        echo "Error: databases are different. See the ipkdiff log"
        exit 6
    fi


    # D140
    D140_REFERENCE="${SCRIPT_DIR}"/data/D140/reference.fasta
    D140_TREE="${SCRIPT_DIR}"/data/D140/tree.newick
    D140_DATABASE_REFERENCE="${SCRIPT_DIR}"/data/D140/DB_k4_o10.ipk
    D140_DATABASE_BUILD="${WORKING_DIR}"/DB_k4_o10.0.ipk

    rm -f "${D140_DATABASE_BUILD}"
    command=python3 "${IPK_SCRIPT}" build -s amino -r "${D140_REFERENCE}" -t "${D140_TREE}" -m LG -k 4 --omega 10.0 -b "${RAXML_NG}" -w "${WORKING_DIR}" -o "${D140_DATABASE_BUILD}"

    echo "Binary files: OK. Running IPK as: ${command}"
    eval "${command}"

    if [ ! -f "${D140_DATABASE_BUILD}" ]
    then
        echo "Error: could not find ${D140_DATABASE_BUILD}. Something went wrong"
        exit 7
    fi

    $IPK_DIFF_AA_BIN 0 "${D140_DATABASE_REFERENCE}" "${D140_DATABASE_BUILD}"

    if [ $? -ne 0 ]; then
        echo "Error: databases are different. See the ipkdiff log"
        exit 8
    fi
fi
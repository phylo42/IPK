#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BIN_DIR="${SCRIPT_DIR}"/../bin
IPK_BIN="${BIN_DIR}"/ipk/ipk-dna
IPK_SCRIPT="${SCRIPT_DIR}"/../ipk.py
IPK_DIFF_BIN="${BIN_DIR}"/tools/ipkdiff-dna
WORKING_DIR="${SCRIPT_DIR}"/../output

REFERENCE="${SCRIPT_DIR}"/data/neotrop/reference.fasta
TREE="${SCRIPT_DIR}"/data/neotrop/tree.rooted.newick
DATABASE_REFERENCE="${SCRIPT_DIR}"/data/neotrop/DB_k7_o2.0.rps
DATABASE_BUILD="${WORKING_DIR}"/DB_k7_o2.0.rps

echo "Script directory: $SCRIPT_DIR"
pushd "${SCRIPT_DIR}"
cd ..

if [ ! -f "${IPK_BIN}" ]
then
    echo "Error: could not find binary files of IPK: ${IPK_BIN}. Please make sure to compile the project"
    exit 1
fi

if [ ! -f "${IPK_DIFF_BIN}" ]
then
    echo "Error: could not find tools: ${IPK_DIFF_BIN}. Please make sure to compile it separately, i.e. do 'make diff-dna' or 'cmake --build DIR --target diff-dna"
    exit 2
fi

if [ ! `which raxml-ng` ]
then
    echo "Error: could not find raxml-ng."
    exit 3
fi

echo "Binary files: OK. Running IPK..."

mkdir -p "${WORKING_DIR}"
rm -f "${DATABASE_BUILD}"
python3 "${IPK_SCRIPT}" build -r "${REFERENCE}" -t "${TREE}" -m GTR -k 7 --omega 2.0 -u 1.0 -b `which raxml-ng` -w "${WORKING_DIR}"

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

popd
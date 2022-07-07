#!/usr/bin/env bash

EXE="${1}"
TESTDIR="${2}"

tfcm1 () {
	# Compare two text files, omitting the first line of each, which is presumed
	# to contain version numbering. A version of the sed editor is required to
	# make a copy of a file, omitting the first line. This function is intended
	# to support tchs.

	# Usage: tfcm1 file1 file2

	# Last revised 04/18/2020 by TJW
	# Converted to bash 2022/07/06 by CRM

	if [ -z "${1}" ]; then
		echo "Error - (${0} tfcm1): First filename is missing. Usage: tfcm1 file1 file2"
		exit 1
	fi

	if [ -z "${2}" ]; then
		echo "Error - (${0} tfcm1): Second filename is missing. Usage: tfcm1 file1 file2"
		exit 1
	fi

	if [ ! -z "${3}" ]; then
		echo "Error - (${0} tfcm1): Third filename is present. Usage: tfcm1 file1 file2"
		exit 1
	fi

	if [ ! -f "${1}" ]; then
		echo "Error - (${0} tfcm1): File "${1}" is not present."
		exit 1
	fi

	if [ ! -f "${2}" ]; then
		echo "Error - (${0} tfcm1): File "${1}" is not present."
		exit 1
	fi

	[ -f acopy1.tmp ] && rm -f acopy1.tmp
	[ -f acopy2.tmp ] && rm -f acopy2.tmp

	# Make a copy of each file, deleting the first line of each.

	sed '1d' "${1}" > acopy1.tmp
	sed '1d' "${2}" > acopy2.tmp

	diff -qZ acopy1.tmp acopy2.tmp > /dev/null 2>&1
	RETVAL="${?}"

	rm -f acopy1.tmp
	rm -f acopy2.tmp

	if [ "${RETVAL}" == 0 ]; then
		echo "The files ${1} and ${2} are identical"
		# set errstr=FALSE
		return 0
	elif [ "${RETVAL}" == 1 ]; then
		echo "The files ${1} and ${2} are different"
		# set errstr=TRUE
		return 1
	fi
}

tchs () {
	# This is a DOS batch file to run an H2OI95 test
	# case. It is used by tch.bat. This batch file in turn
	# uses tfcm1.bat to compare output files. This particular
	# comparison disregards the first line, which contains
	# the H2OI95 build information.

	# Last revised 05/18/2020 by TJW
	# Converted to bash 2022/07/06 by CRM

	# First argument is the absolute location of the executable
	EXE="${1}"

	# Here "${2}" is the name of the specified test case
	# folder, which corresponds a single test case.
	TEST=${2}

	if [ -z "${TEST}" ]; then
		echo "Error - (${0} tchs): Missing argument."
		exit 1
	fi

	if [ ! -d "${TEST}" ]; then
		echo "Error - (${0} tchs): Folder ${TEST} is not present."
		exit 1
	fi

	# Change to the test case folder.
	cd "${TEST}"
	
	TEST=`basename ${TEST}`
	echo "Processing test case ${TEST}"

	# A valid test case folder must contain a file
	# named input, or an archival copy thereof.

	# Get the test case string without the leading
	# two characters, which denote the test case type.
	STR="${TEST:2}"
	echo "The short string is ${STR}"

	INPUTA="input_${STR}"

	echo "Copying ${INPUTA} to input"
	cp "${INPUTA}" input

	OUTPUTA="output_${STR}"
	OUTA="out_${STR}"
	MTABA="mtab_${STR}.csv"
	CTABA="ctab_${STR}.csv"
	XTABA="xtab_${STR}.csv"

	echo "Running H2OI95"
	"${EXE}" > out

	DIFFS=false

	# Compare the normal and archival output files

	tfcm1 output "${OUTPUTA}" || DIFFS=true

	# Compare the normal and archival out files

	tfcm1 out "${OUTA}" || DIFFS=true

	# Compare the normal and archival mtab files

	tfcm1 mtab.csv "${MTABA}" || DIFFS=true

	# Compare the normal and archival ctab files

	tfcm1 ctab.csv "${CTABA}" || DIFFS=true

	# Compare the normal and archival xtab files

	tfcm1 xtab.csv "${XTABA}" || DIFFS=true

	UPDATE_ARCHIVAL=false
	REMOVE_OUTPUT=true

	if [ "${UPDATE_ARCHIVAL}" == true ];then
		# The following if active updates archival output files.
		# This is not the normal testing procedure.

		echo "Update archival output files"

		mv -f output "${OUTPUTA}"
		mv -f out "${OUTA}"
		mv -f mtab.csv "${MTABA}"
		mv -f ctab.csv "${CTABA}"
		mv -f xtab.csv "${XTABA}"
	fi
	if [ "${REMOVE_OUTPUT}" == true ]; then
		# The following line if active skips removing the normal input and output files.
		# This is the normal testing procedure (removing).

		rm -f input output out mtab.csv ctab.csv xtab.csv
	fi

	# Return to the parent folder.
	cd - > /dev/null 2>&1

	if [ "${DIFFS}" == true ]; then
		echo "----- DIFFERENCES DETECTED -----"
		if [ "${UPDATE_ARCHIVAL}" == false ]; then
			return 1
		else
			echo "----- OUTPUTS UPDATED -----"
			return 0
		fi
	else
		echo "NO DIFFERENCES FOUND"
		return 0
	fi
}

tchs "${EXE}" "${TESTDIR}"
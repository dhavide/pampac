# Modify these variables as appropriate for your system.
MPIEXEC := /usr/bin/mpiexec
# Default value of NPROC is 4; can be set at command line
NPROC ?= 4
DATADIR := example/data
PARAMETERFILE := example/parameters.txt
INPUTFILE := ${DATADIR}/input.txt_2048
OUTPUTFILE := ${DATADIR}/output.txt

#
TARGET := ${OUTPUTFILE}
TGT_PREREQS := example/KS_example
TGT_POSTMAKE := 
TGT_POSTCLEAN := rm -f ${OUTPUTFILE} ${DATADIR}/*.gv ${DATADIR}/*.svg
GV_FILES := $(wildcard ${DATADIR}/*.gv)
SVG_FILES := $(GV_FILES:.gv=.svg)

KS_exec: 
	${MPIEXEC} -n ${NPROC} ${TGT_PREREQS} ${PARAMETERFILE} ${INPUTFILE} ${OUTPUTFILE} ${DATADIR}/tree

%.svg : %.gv
	dot -Tsvg $< -o $*.svg

svg: $(SVG_FILES)

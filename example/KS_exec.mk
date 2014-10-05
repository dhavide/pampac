# Modify these variables as appropriate for your system.
MPIEXEC := /usr/bin/mpiexec
# Default value of NPROC is 4; can be set at command line
NPROC ?= 4
DATADIR := example/data
PARAMETERFILE := example/parameters.txt

TGT_PREREQS := example/KS_example
TGT_POSTMAKE := 
TGT_POSTCLEAN := rm -f ${DATADIR}/*.gv ${DATADIR}/*.svg
GV_FILES := $(wildcard ${DATADIR}/*.gv)
SVG_FILES := $(GV_FILES:.gv=.svg)

KS_exec: 
	${MPIEXEC} -n ${NPROC} ${TGT_PREREQS} ${PARAMETERFILE} 

%.svg : %.gv
	dot -Tsvg $< -o $*.svg

svg: $(SVG_FILES)

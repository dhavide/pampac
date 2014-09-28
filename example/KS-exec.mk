TARGET := output.txt

TGT_PREREQS := KS-example
TGT_POSTMAKE := 
TGT_POSTCLEAN := rm -f example/output.txt example/tree-data/*.gv example/tree-data/*.svg
MPIEXEC := /usr/bin/mpiexec
NPROC ?= 4

GV_FILES := $(wildcard example/tree-data/*.gv)
SVG_FILES := $(GV_FILES:.gv=.svg)

KS-exec: example/KS-example
	${MPIEXEC} -n ${NPROC} example/KS-example example/parameters.txt \	
	    example/input.txt example/output.txt example/tree-data/tree

%.svg : %.gv
	dot -Tsvg $< -o $*.svg

svg: $(SVG_FILES)

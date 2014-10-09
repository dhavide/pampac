
.PHONY: trees
TGT_POSTCLEAN := rm -f ${DATADIR}/*.gv ${DATADIR}/*.svg
GV_FILES := $(wildcard ${DATADIR}/*.gv)
SVG_FILES := $(GV_FILES:.gv=.svg)


%.svg : %.gv
	dot -Tsvg $< -o $*.svg

trees: $(SVG_FILES)

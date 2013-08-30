GV_FILES := $(wildcard *.gv)
SVG_FILES := $(notdir $(GV_FILES:.gv=.svg))

.SUFFIXES:
.SUFFIXES: .gv .svg

svg: $(SVG_FILES)

%.svg: %.gv
	dot -Tsvg $< -o $*.svg

clean:
	@rm -f $(SVG_FILES) $(GV_FILES)

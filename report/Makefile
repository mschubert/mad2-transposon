RSCRIPTS = $(wildcard *.r)
SVG = $(wildcard external/*.svg)
PDF = $(RSCRIPTS:%.r=%.pdf)
PNG = $(PDF:%.pdf=%.png)
SVGPNG = $(SVG:%.svg=%.png)
PPSVG = $(wildcard postproc/*.svg)
PPPDF = $(PPSVG:%.svg=%.pdf)

all: $(PDF) $(PNG) $(PPPDF)

%.pdf: %.r
	Rscript $^ $@

%.png: %.pdf
	convert -flatten -density 400 $< -resize 25% $@

$(PPPDF): %.pdf: %.svg
	inkscape $^ --export-pdf=$@

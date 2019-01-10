RSCRIPTS = $(wildcard *.r)
SVG = $(wildcard external/*.svg)
PDF = $(RSCRIPTS:%.r=%.pdf)
PNG = $(PDF:%.pdf=%.png)
SVGPNG = $(SVG:%.svg=%.png)

all: $(PDF) $(PNG)

%.pdf: %.r
	Rscript $^ $@

%.png: %.pdf
	convert -flatten -density 400 $< -resize 25% $@

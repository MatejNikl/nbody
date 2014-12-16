GRAPH_DATA := $(wildcard data/*)
GP_SRC     := $(sort $(wildcard *.gp))
GRAPHS     := $(GP_SRC:.gp=.eps)

TEX_SRC    := $(sort $(wildcard *.tex))
PDFS       := $(TEX_SRC:.tex=.pdf)


all: graph pdf

graph: $(GRAPHS)

pdf: $(PDFS)

%.eps: %.gp $(GRAPH_DATA)
	gnuplot $<

%.pdf: %.tex
	pdflatex $^

.PHONY: all graph pdf

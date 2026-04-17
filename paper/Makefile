# Makefile for Paper_v1 (variance-corrected SIM composition)
# Requires: latexmk, pdflatex, bibtex

PAPER   := Paper_v1
BIB     := References_v1
TEX     := $(PAPER).tex
PDF     := $(PAPER).pdf
SECTIONS := $(wildcard sections/*.tex) $(wildcard sections/tables/*.tex)
FIGS     := $(wildcard figs/*)

LATEXMK := latexmk
LATEXMK_FLAGS := -pdf -bibtex -interaction=nonstopmode -halt-on-error

.PHONY: all pdf watch clean distclean view

all: pdf

pdf: $(PDF)

$(PDF): $(TEX) $(BIB).bib $(SECTIONS) $(FIGS) neurips_2026.sty
	$(LATEXMK) $(LATEXMK_FLAGS) $(PAPER)

# Continuous rebuild on file changes (Ctrl-C to stop).
watch:
	$(LATEXMK) $(LATEXMK_FLAGS) -pvc $(PAPER)

# Open the built PDF (macOS).
view: $(PDF)
	open $(PDF)

# Remove intermediate build artifacts, keep the PDF.
clean:
	$(LATEXMK) -c $(PAPER)
	rm -f $(PAPER).bbl $(PAPER).run.xml $(PAPER).synctex.gz

# Remove everything, including the PDF.
distclean:
	$(LATEXMK) -C $(PAPER)
	rm -f $(PAPER).bbl $(PAPER).run.xml $(PAPER).synctex.gz

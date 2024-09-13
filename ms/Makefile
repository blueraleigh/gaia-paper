.PHONY: gaia

gaia: gaia.pdf

gaia.pdf: gaia.tex gaia.bib
	pdflatex gaia.tex
	bibtex gaia.aux
	pdflatex gaia.tex
	pdflatex gaia.tex

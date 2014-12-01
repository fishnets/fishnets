all: steepness_empirical.pdf

steepness_empirical.pdf: steepness_empirical.tex
	xelatex steepness_empirical.tex; biber steepness_empirical; xelatex steepness_empirical.tex

steepness_empirical.tex: steepness_empirical.Rnw
	Rscript --vanilla -e 'require(knitr); knit("steepness_empirical.Rnw")'

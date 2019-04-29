#!/bin/bash

# change into requested directory
cd $1

# clean
rm -vf main.aux
rm -vf main.bbl
rm -vf main.bcf
rm -vf main.blg
rm -vf main.log
rm -vf main.out
rm -vf main.pdf
rm -vf main.run.xml
rm -vf main.synctex.gz

# ensure bail out on error
set -eo pipefail

# run pdflatex + biber + pdflatex
pdflatex -interaction=batchmode -halt-on-error main.tex
biber main
pdflatex -interaction=batchmode -halt-on-error main.tex

# return to parent directory
cd ..

#!/bin/bash

# ensure bail out on error
set -eo pipefail

# setup environment
source binder/env.sh

function rebuild {

    # clean
    rm -vf ${1}.aux
    rm -vf ${1}.bbl
    rm -vf ${1}.bcf
    rm -vf ${1}.blg
    rm -vf ${1}.log
    rm -vf ${1}.out
    rm -vf ${1}.lof
    rm -vf ${1}.lot
    rm -vf ${1}.toc
    rm -vf ${1}.pdf
    rm -vf ${1}.run.xml
    rm -vf ${1}.synctex.gz

    # run pdflatex + biber + pdflatex
    pdflatex -interaction=nonstopmode -halt-on-error ${1}.tex
    biber ${1}
    pdflatex -interaction=nonstopmode -halt-on-error ${1}.tex

}

rebuild abstract
rebuild chapter1

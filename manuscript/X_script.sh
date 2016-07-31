#/bin/bash
LATEXFILE=$1

filename=$(basename "$LATEXFILE")  # from http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
extension="${filename##*.}"
filename="${filename%.*}"

latex $LATEXFILE
bibtex $filename
latex $LATEXFILE
latex $LATEXFILE
dvips $filename
ps2pdf $filename.ps

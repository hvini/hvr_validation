#!/usr/bin/env bash

## Shell script to generate a single LaTeX file compiling all tables
## and compile it

# HERE
HERE=/ipl/ipl27/sfernandez/hvr_validation
TABLESDIR=${HERE}/tables
OUTPUT=${HERE}/tables/tables.tex

## Create the master file with the preamble
cat << 'EOF' > $OUTPUT
\documentclass{article}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{multirow}
\begin{document}
EOF

# Loop through all table files (adjust the pattern as needed)
for file in ${TABLESDIR}/adni-bl*.tex
do
	#echo "\section*{$title}" >> $OUTPUT
	echo "\input{$file}" >> $OUTPUT
	echo "\clearpage" >> $OUTPUT

	done

	# Close the document
	echo "\end{document}" >> $OUTPUT

	# Optionally compile the document
	pdflatex -output-directory=$TABLESDIR $OUTPUT

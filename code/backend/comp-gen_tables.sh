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
for file in ${TABLESDIR}/adni-bl_table-?.tex
do
	# Remove stub line
	sed -i -e '/l|c/ s/|//' $file
	# Add title
	if [[ $file =~ ([0-9])\.tex$ ]]
	then
		tablenum="${BASH_REMATCH[1]}"
		case "$tablenum" in
			1) title="Table 1: Demographic data" ;;
			2) title="Table 2: HCvol, HVR, and segmentation failures" ;;
			*) title="Table $tablenum" ;;
		esac
		echo "\section*{$title}" >> $OUTPUT
	fi

	# Insert table
	echo "\input{$file}" >> $OUTPUT

	# Flush
	echo "\clearpage" >> $OUTPUT

	done

	# Close the document
	echo "\end{document}" >> $OUTPUT

	# Optionally compile the document
	pdflatex -output-directory=$TABLESDIR $OUTPUT

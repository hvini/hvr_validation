#!/usr/bin/env bash

## Shell script to generate a single LaTeX file compiling all tables

## CONSTANTS
HERE=/ipl/ipl27/sfernandez/hvr_validation
TABLESDIR=${HERE}/tables
METADATA=${TABLESDIR}/metadata.json
OUTPUT=${TABLESDIR}/tables.tex

## Create the master file with the preamble
cat << 'EOF' > $OUTPUT
\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{longtable}
\usepackage{rotating}
\usepackage{multirow}
\usepackage{caption}
\usepackage[flushleft,para]{threeparttable}
\begin{document}
\begin{titlepage}
\begin{center}
\vspace*{3cm}\par
\Large
\textbf{Enhanced Detection of Age-Related and Cognitive Declines Using Automated Hippocampal-to-Ventricle Ratio in Alzheimer's Patients}
\vspace*{2.5cm}\par
\Huge
Tables
\end{center}
\end{titlepage}
\thispagestyle{empty}
\listoftables
\clearpage
\pagenumbering{arabic}
\captionsetup{width=0.8\textwidth,justification=centering,labelfont={large,sc},labelsep=endash,textfont=normalsize}
EOF

# Loop through all table files (adjust the pattern as needed)
for file in ${TABLESDIR}/table-[0-9]*.tex
do
	## Capture table number for purposes...
	if [[ $file =~ table-([0-9]) ]]
	then
		tablenum="${BASH_REMATCH[1]}"
	else
		continue
	fi

	# Name of modified file that will be used as input
	outfile=${file/table-/table_clean-}

	# Remove empty (last) line
	sed '/^$/d' $file > $outfile

	# Go from table to threeparttable
	sed -i \
		-e '1s/.*/\\begin{threeparttable}/' \
		-e '$s/.*/\\end{threeparttable}/' \
		$outfile

	# Remove stub line & change tabular* to tabular
	sed -i -E '3s/.*l\|(c+).*/\\begin{tabular}{l\1}/' $outfile
	sed -i -E 's/\\end\{tabular\*\}/\\end{tabular}/' $outfile

	# Deal with minipage/footnote nonsense
	sed -i -E \
		-e 's/\\begin\{minipage.*/\\begin{tablenotes}\\small/' \
		-e 's/\\end\{minipage.*/\\end{tablenotes}/' \
		$outfile

	# Caption
	# Use jq to extract the title and caption from metadata.json
	title=$(jq -r ".table${tablenum}.title" $METADATA)
	caption=$(jq -r ".table${tablenum}.caption" $METADATA)

	sed -i "/\\end{threeparttable}/i\\
\\\\caption[$title]{\\
	\\\\textbf{$title}\\\\\\\\$caption\\\\\\\\}" $outfile

	# Insert table
	echo "\centering" >> $OUTPUT
	echo "\input{$outfile}" >> $OUTPUT
	if [[ "$tablenum" -eq 1 ]]
	then
		echo "\clearpage" >> $OUTPUT # Flush new page
	fi
done

# Close the document
echo "\end{document}" >> $OUTPUT

# Optionally compile the document
pdflatex -output-directory=$TABLESDIR $OUTPUT
pdflatex -output-directory=$TABLESDIR $OUTPUT

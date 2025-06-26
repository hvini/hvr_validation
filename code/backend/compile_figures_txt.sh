#!/usr/bin/env bash

## Shell script to generate a single LaTeX file compiling all figures

## CONSTANTS
HERE=/ipl/ipl27/sfernandez/hvr_validation
FIGSDIR=${HERE}/plots
METADATA=${FIGSDIR}/metadata.json
OUTPUT=${FIGSDIR}/figures.tex

## Create the master file with the preamble
cat << 'EOF' > $OUTPUT
\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{caption}
\begin{document}
\begin{titlepage}
\begin{center}
\vspace*{3cm}\par
\Large
\textbf{Enhanced Detection of Age-Related and Cognitive Declines Using Automated Hippocampal-to-Ventricle Ratio in Alzheimer's Patients}
\vspace*{2.5cm}\par
\Huge
Figures
\end{center}
\end{titlepage}
\thispagestyle{empty}
\listoffigures
\clearpage
\pagenumbering{arabic}
\captionsetup{width=0.8\textwidth,justification=centering,labelfont={large,sc},labelsep=endash,textfont=normalsize}
EOF

# Loop through all table files (adjust the pattern as needed)
for file in ${FIGSDIR}/fig-[0-9]_*.png
do
	## Capture table number for purposes...
	if [[ $file =~ fig-([0-9]) ]]
	then
		fignum="${BASH_REMATCH[1]}"
	else
		continue
	fi

	# Caption
	# Use jq to extract the title and caption from metadata.json
	title=$(jq -r ".fig${fignum}.title" $METADATA)
	caption=$(jq -r ".fig${fignum}.caption" $METADATA)

	# Insert figure
	echo "\begin{figure}" >> $OUTPUT
	echo "\centering" >> $OUTPUT
	if [[ "$fignum" -eq 2 ]]
	then
		echo "\includegraphics{$file}" >> $OUTPUT
	else
		echo "\includegraphics[width=\textwidth]{$file}" >> $OUTPUT
	fi
	echo "\caption[$title]{\textbf{$title}\\\\$caption\\\\}" >> $OUTPUT
	echo "\end{figure}" >> $OUTPUT
	echo "\clearpage" >> $OUTPUT
done

# Close the document
echo "\end{document}" >> $OUTPUT

# Optionally compile the document
pdflatex -output-directory=$FIGSDIR $OUTPUT
pdflatex -output-directory=$FIGSDIR $OUTPUT

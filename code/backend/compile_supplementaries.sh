#!/usr/bin/env bash

## Shell script to generate a single LaTeX file compiling all tables
## and compile it

# HERE
HERE=/ipl/ipl27/sfernandez/hvr_validation
TABLESDIR=${HERE}/tables
FIGSDIR=${HERE}/plots
TMETADATA=${TABLESDIR}/metadata.json
FMETADATA=${FIGSDIR}/metadata.json
OUTDIR=${HERE}/supplementaries
OUTPUT=${OUTDIR}/supplementaries.tex

## Create output directory
[[ -d $OUTDIR ]] || mkdir $OUTDIR

## Create the master file with the preamble
cat << 'EOF' > $OUTPUT
\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{longtable}
\usepackage{rotating}
\usepackage{multirow}
\usepackage[flushleft,para]{threeparttable}
\usepackage{graphicx}
\usepackage{caption}
\DeclareCaptionLabelFormat{sup}{#1 S#2}
\captionsetup{width=0.8\textwidth,justification=centering,labelformat=sup,labelfont={large,sc},labelsep=endash,textfont=normalsize}
\begin{document}
\begin{titlepage}
\begin{center}
\vspace*{3cm}\par
\Large
\textbf{Enhanced Detection of Age-Related and Cognitive Declines Using Automated Hippocampal-to-Ventricle Ratio in Alzheimer's Patients}
\vspace*{2cm}\par
\Huge
Supplementary material
\end{center}
\end{titlepage}
\thispagestyle{empty}
\renewcommand{\listtablename}{List of Supplementary Tables}
\renewcommand{\listfigurename}{List of Supplementary Figures}
\listoftables
\listoffigures
\clearpage
\pagenumbering{arabic}
EOF

### Tables
echo "\centering" >> $OUTPUT
echo "\section*{Supplementary Tables}" >> $OUTPUT
echo "\vspace*{2cm}\par" >> $OUTPUT
#echo "\captionsetup[table]{name=Supplementary Table}" >> $OUTPUT
# Loop through all table files (adjust the pattern as needed)
for file in ${TABLESDIR}/table-s[0-9]*.tex
do
	## Capture table number for purposes...
	if [[ $file =~ table-s([0-9]) ]]
	then
		tablenum="${BASH_REMATCH[1]}"
	else
		continue
	fi

	# Name of modified file that will be used as input
	outfile=${file/table-s/table_clean-s}

	# Remove empty lines
	sed '/^$/d' $file > $outfile

	# Go from table to threeparttable
	sed -i \
		-e '1s/.*/\\begin{threeparttable}/' \
		-e '$s/.*/\\end{threeparttable}/' \
		$outfile

	# Fontsize
	case "$tablenum" in
		4|5) sed -i -E '2s/.*/\\fontsize{10.0pt}{12.0pt}\\selectfont/' "$outfile";;
	esac

	## Remove stub line
	case "$tablenum" in
		4|5) sed -i -E '3s/(l)\|(c+)/c\2/' "$outfile";;
		*) sed -i -E '3s/(l)\|(c+)/\1\2/' "$outfile";;
	esac

	## Remove stub line & change tabular* to tabular
	#sed -i -E '3s/.*l\|(c+).*/\\begin{tabular}{l\1}/' $outfile
	#sed -i -E 's/\\end\{tabular\*\}/\\end{tabular}/' $outfile

	# Deal with minipage/footnote nonsense
	sed -i -E \
		-e 's/\\begin\{minipage.*/\\begin{tablenotes}\\small/' \
		-e 's/\\end\{minipage.*/\\end{tablenotes}/' \
		$outfile

	# Caption
	# Use jq to extract the title and caption from metadata.json
	title=$(jq -r ".table${tablenum}s.title" $TMETADATA)
	caption=$(jq -r ".table${tablenum}s.caption" $TMETADATA)

	sed -i "/\\end{threeparttable}/i\\
\\\\caption[$title]{\\
	\\\\textbf{$title}\\\\\\\\$caption\\\\\\\\}" $outfile

	# Insert table
	echo "\centering" >> $OUTPUT
	echo "\input{$outfile}" >> $OUTPUT
	#
	# Flush
	#case "$tablenum" in
		#1|4) echo "\clearpage" >> $OUTPUT;;
	#esac
done

### Figures
echo "\clearpage" >> $OUTPUT
echo "\section*{Supplementary Figures}" >> $OUTPUT
echo "\vspace*{2cm}\par" >> $OUTPUT
#echo "\captionsetup[figure]{name=Supplementary Figure}" >> $OUTPUT

# Loop through all table files (adjust the pattern as needed)
for file in ${FIGSDIR}/fig-s[0-9]*.png
do
	## Capture table number for purposes...
	if [[ $file =~ fig-s([0-9]) ]]
	then
		fignum="${BASH_REMATCH[1]}"
	else
		continue
	fi

	# Caption
	# Use jq to extract the title and caption from metadata.json
	title=$(jq -r ".fig${fignum}s.title" $FMETADATA)
	caption=$(jq -r ".fig${fignum}s.caption" $FMETADATA)

	# Insert figure
	echo "\begin{figure}[h]" >> $OUTPUT
	echo "\centering" >> $OUTPUT
	## TODO: Check which require size adjustment
	case "$fignum" in
		1|2|3) echo "\includegraphics{$file}" >> "$OUTPUT";;
		*) echo "\includegraphics[width=\textwidth]{$file}" >> "$OUTPUT";;
	esac
	echo "\caption[$title]{\textbf{$title}\\\\$caption\\\\}" >> $OUTPUT
	echo "\end{figure}" >> $OUTPUT
	echo "\clearpage" >> $OUTPUT
done

# Close the document
echo "\end{document}" >> $OUTPUT

# Optionally compile the document
pdflatex -output-directory=$OUTDIR $OUTPUT
pdflatex -output-directory=$OUTDIR $OUTPUT

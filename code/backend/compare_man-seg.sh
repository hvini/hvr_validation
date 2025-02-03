#!/usr/bin/env bash

## Compare overlap similarity between segmentations:
## { MALF, NLPB, CNN } — manual labels
## XCorrelation

TMPDIR=$(mktemp -d --tmpdir)
trap "rm -rf $TMPDIR" 0 1 2 15

set -ux

## FUNCTIONS
# Recode labels to L-HC: 1 & R-HC: 2 & L-VC: 3 & R-VC: 4
# Reshape to Dorothee's labels
recode() {
	local input="$1"
	local reference="$2"
	local output="$3"
	itk_resample \
		"$input" \
		--labels \
		--clobber \
		--lut-string "11 1; 12 3; 21 2; 22 4" \
		--like "$reference" \
		"$output"
}

# Compare two label volumes
compare() {
	local subj_id="$1"
	local img1="$2"
	local img2="$3"
	minccmp -similarity "$img1" "$img2" | awk -v id="$subj_id" '
	$1 ~ /^[1234]$/ {
		# Left/Right
		if ($1 == "1" || $1 == "3")
			side = "left";
		else if ($1 == "2" || $1 == "4")
			side = "right";
		# HC/VC
		if ($1 == "1" || $1 == "2")
			roi = "hc";
		else if ($1 == "3" || $1 == "4")
			roi = "vc";
		dice = $2
		sens = $3
		spec = $4
		accu = $5
		kppa = $6
		print id "," roi "," side "," dice "," kppa "," accu "," sens "," spec
	}'
}

# PATHS
HERE=/ipl/ipl27/sfernandez/hvr_validation
MAN_LABS=${HERE}/data/labels_dorothee/labels_reduced
MALF_LABS=${HERE}/data/orig_segmentations/malf/reduced
NLPB_LABS=${HERE}/data/orig_segmentations/nlpb/reduced
CNN_LABS=${HERE}/data/orig_segmentations/cnn/hcvc
OUTFILE=${HERE}/data/derivatives/man-seg_kappa_hcvc_method.csv

CSV_HEADER="id,roi,side,dice,kappa,accuracy,sensitivity,specificity"
for segm in malf nlpb cnn
do
	echo $CSV_HEADER > ${OUTFILE/method/$method}
done

for manual_lab in ${MAN_LABS}/*
do
	# Extract id
	id=$(echo $manual_lab | grep -oP "0\d{2}")

	# Reshape manual label
	man_=${TMPDIR}/${id}_man_recoded.mnc
	recode "$manual_lab" "$manual_lab" "$man_"

	# MALF
	# Find corresponding segmented labels
	malf=$(find $MALF_LABS -name "*$id*")

	# Recode & reshape
	malf_=${TMPDIR}/${id}_malf_reshaped.mnc
	recode "$malf"  "$manual_lab" "$malf_"

	# Comparison
	compare $id $man_ $malf_ >> ${OUTFILE/method/malf}

	# NLPB
	# Find corresponding segmented labels
	nlpb=$(find $NLPB_LABS -name "*$id*")

	# Recode & reshape
	nlpb_=${TMPDIR}/${id}_nlpb_reshaped.mnc
	recode "$nlpb"  "$manual_lab" "$nlpb_"

	# Comparison
	compare $id $man_ $nlpb_ >> ${OUTFILE/method/nlpb}

	# CNN
	# Find corresponding segmented labels
	cnn=$(find $CNN_LABS -name "*$id*")

	# Recode & reshape
	cnn_=${TMPDIR}/${id}_cnn_reshaped.mnc
	recode "$cnn"  "$manual_lab" "$cnn_"

	# Comparison
	compare $id $man_ $cnn_ >> ${OUTFILE/method/cnn}
done

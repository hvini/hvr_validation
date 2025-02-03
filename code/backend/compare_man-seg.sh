#!/usr/bin/env bash

## Compare overlap similarity between segmentations:
## { MALF, NLPB, CNN } — manual labels
## XCorrelation

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

# Process
label_comparison() {
	local man_lab="$1"
	local TMPDIR=$(mktemp -d --tmpdir)
	trap "rm -rf $TMPDIR" EXIT

	# Extract id
	local id=$(echo $man_lab | grep -oP "0\d{2}")

	# Reshape manual label
	local man_lab_=${TMPDIR}/${id}_man_recoded.mnc
	recode "$man_lab" "$man_lab" "$man_lab_"

	# Process each segmentation method
	for segm in malf nlpb cnn
	do
		local path_var="${segm^^}_DIR"
		local seg_lab=$(find "${!path_var}" -name "*${id}*")
		local seg_lab_="${TMPDIR}/${id}_${segm}_reshaped.mnc"
		recode "$seg_lab" "$man_lab_" "$seg_lab_"
		compare "$id" "$man_lab_" "$seg_lab_" >> "${OUTFILE/method/$segm}"
	done
}

# Export functions
export -f recode
export -f compare
export -f label_comparison

# CONSTANTS
HERE=/ipl/ipl27/sfernandez/hvr_validation
MAN_DIR=${HERE}/data/labels_dorothee/labels_reduced
MALF_DIR=${HERE}/data/orig_segmentations/malf/reduced
NLPB_DIR=${HERE}/data/orig_segmentations/nlpb/reduced
CNN_DIR=${HERE}/data/orig_segmentations/cnn/hcvc
OUTFILE=${HERE}/data/derivatives/man-seg_kappa_hcvc_method.csv

# Export constants
export MALF_DIR NLPB_DIR CNN_DIR OUTFILE

# INIT csv files
CSV_HEADER="id,roi,side,dice,kappa,accuracy,sensitivity,specificity"
for segm in malf nlpb cnn
do
	echo $CSV_HEADER > ${OUTFILE/method/$segm}
done

## Parallel comparison
find "${MAN_DIR}" -type f | parallel --jobs 6 label_comparison

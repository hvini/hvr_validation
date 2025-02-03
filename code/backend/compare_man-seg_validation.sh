#!/usr/bin/env bash

## Compare overlap similarity between segmentations:
## CNN and manual labels :: Validation datasets — ADNI & ICBM
## XCorrelation

set -ux

## HOME
HERE=/ipl/ipl27/sfernandez/hvr_validation

## FUNCTIONS
# Recode labels to L-HC: 1 & R-HC 2
recode() {
	local input="$1"
	local reference="$2"
	local output="$3"
	local left_hc="$4"
	local right_hc="$5"
	printf -v lut_string "%i 1; %i 2" "$left_hc" "$right_hc"
	itk_resample \
		"$input" \
		--labels \
		--clobber \
		--lut-string "$lut_string" \
		--like "$reference" \
		"$output"
}

# Compare two label volumes
compare() {
	local subj_id="$1"
	local img1="$2"
	local img2="$3"
	minccmp -similarity "$img1" "$img2" | awk -v id="$subj_id" '
	$1 ~ /^[12]$/ {
		if ($1 == "1")
			hc = "hc_left";
		else if ($1 == "2")
			hc = "hc_right";
		dice = $2
		sens = $3
		spec = $4
		accu = $5
		kppa = $6
		print id "," hc "," dice "," kppa "," accu "," sens "," spec
	}'
}


## ICBM
# CONSTANTS
ICBM=${HERE}/data/validation/jens_icbm
ICBM_OUTDIR=${ICBM}/comparison_cnn_manual
ICBM_OUTFILE=${ICBM_OUTDIR}/comparison_cnn_manual_icbm.csv

[ -d "$ICBM_OUTDIR" ] || mkdir $ICBM_OUTDIR

ICBM_MANUAL=${ICBM}/labels
ICBM_CNN=${ICBM}/segm/segmentations

printf "id,hc,dice,kappa,accuracy,sensitivity,specificity\n" > $ICBM_OUTFILE

for man_lab in ${ICBM_MANUAL}/*
do
	# Extract subject_id
	id=$(echo $man_lab | grep -oP "\d{5}" -)

	# Recode manual label
	out_man=${ICBM_OUTDIR}/hc_icbm_man_${id}.mnc
	# Labels: LHC 4 & RHC 2
	recode "$man_lab" "$man_lab" "$out_man" 4 2

	# Recode & reshape cnn label
	cnn_lab=$(find $ICBM_CNN -name "*${id}*")
	out_cnn=${out_man/_man_/_cnn_}
	# Labels: LHC 11 & RHC 21
	recode "$cnn_lab" "$out_man" "$out_cnn" 11 21

	# Comparison
	compare $id $out_man $out_cnn >> $ICBM_OUTFILE
done


## ADNI
ADNI=${HERE}/data/validation/jens_adni
ADNI_OUTDIR=${ADNI}/comparison_cnn_manual
ADNI_OUTFILE=${ADNI_OUTDIR}/comparison_cnn_manual_adni.csv

[ -d "$ADNI_OUTDIR" ] || mkdir $ADNI_OUTDIR

ADNI_MANUAL=${ADNI}/labels_manual_src
ADNI_CNN=${ADNI}/segm/segmentations

printf "id,hc,dice,kappa,accuracy,sensitivity,specificity\n" > $ADNI_OUTFILE

# Because I need to match IDs, this one uses the CNN labels
for cnn_lab in ${ADNI_CNN}/*
do
	# Extract subject_ids (original & labels')
	id1=$(echo $cnn_lab | grep -oP "\d+_S_\d+")
	id2=$(minchistory $cnn_lab | grep -oPm 1 "ADNI\d{3}")

	# Find corresponding manual label
	man_lab=$(find $ADNI_MANUAL -name "*${id2:(-3)}*")

	# Recode manual label
	out_man=${ADNI_OUTDIR}/hc_adni_man_${id1}.mnc
	# Labels: LHC 16 & RHC 3
	recode "$man_lab" "$man_lab" "$out_man" 16 3

	# Recode & reshape cnn label
	out_cnn=${out_man/_man_/_cnn_}
	# Labels: LHC 11 & RHC 21
	recode "$cnn_lab" "$out_man" "$out_cnn" 11 21

	# Comparison
	compare $id1 $out_man $out_cnn >> $ADNI_OUTFILE
done

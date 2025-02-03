#!/usr/bin/env bash

### Parse HCVC segmentations (simple)
## Extract voxel num for each label and calculate complete HV/CSF
## Labels:
## 11 - L-HC
## 12 - L-CSF
## 21 - R-HC
## 22 - R-CSF

set -ex

IN_DIR=$1
OUT_FILE=$2
LHC=${3:-11}
RHC=${4:-21}
LCSF=${5:-12}
RCSF=${6:-22}

[ $# -ne 2 ] || [ ! -d $IN_DIR ] && exit 1

[ -f $OUT_FILE ] && rm $OUT_FILE
printf "ID,LHC,RHC,HC,LCSF,RCSF,CSF\n" > $OUT_FILE

for img in $IN_DIR/*
do
	bname=$(basename $img .mnc)

	lhc=$(print_all_labels $img | awk -v lhc=$LHC '/Label: lhc/ {print $3}')
	rhc=$(print_all_labels $img | awk -v rhc=$RHC '/Label: rhc/ {print $3}')
	[ -z $lhc ] && lhc=0
	[ -z $rhc ] && rhc=0
	hc=$(($lhc + $rhc))

	lcsf=$(print_all_labels $img | awk -v lcsf=$LCSF '/Label: lcsf/ {print $3}')
	rcsf=$(print_all_labels $img | awk -v rcsf=$RCSF '/Label: rcsf/ {print $3}')
	[ -z $csf ] && lcsf=0
	[ -z $rcsf ] && rcsf=0
	csf=$(($lcsf + $rcsf))
	printf "%s,%i,%i,%i,%i,%i,%i\n" \
		$bname $lhc $rhc $hc $lcsf $rcsf $csf >> $OUT_FILE
done

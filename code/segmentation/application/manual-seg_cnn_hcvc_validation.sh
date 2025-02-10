#! /usr/bin/env bash

## Apply CNN ensemble models to subjects from validation datasets
## Need to load hvr_validation environment

HERE=/ipl/ipl27/sfernandez/hvr_validation

for dataset in jens_adni jens_icbm
do
datadir=${HERE}/data/validation/${dataset}/t1
outdir=${HERE}/data/validation/${dataset}/cnn
[ -d $outdir ] || mkdir -p $outdir
inputfile=${datadir}/input.csv
if [ -f $inputfile ]
then
	echo printf "Error: File %s not found.\n" $inputfile >&2
	exit 1
fi

docker \
	run --rm \
	-v ${datadir}:/app/data \
	-v ${outdir}:/app/output \
	--user $(id -u) \
	hvr_cnn:latest -f input.csv --has_header --vols --qc

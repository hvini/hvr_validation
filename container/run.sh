#!/bin/bash

docker run --rm --gpus all \
	-v $(pwd)/data:/app/data \
	-v $(pwd)/output:/app/output \
	vhpcavalcanti/hvr_cnn_gpu:latest \
	-i /app/data/IXI002-Guys-0828-T1.mnc \
	-o /app/output --qc

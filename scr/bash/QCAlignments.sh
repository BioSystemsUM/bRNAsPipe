#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline
# It produces graphs for alignments' evaluation

### Variables:
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='quay.io/biocontainers/qualimap:2.2.2d--1'
BaseDir=$1
AlignmentsDir=$2
QCAlignments=$3
Annotation=$4
### Create Qualimap report files:

DockerId=$($docker run --rm --privileged -d -e AlignmentsDir=$AlignmentsDir -e QCAlignments=$QCAlignments -e Annotation=$Annotation -v $BaseDir:/data $Image \
bash -c 'cd /data/$AlignmentsDir ; \
for file in $( ls | grep -E ".*Sorted.bam$")
do DirName=$(echo $file | sed "s/.Sorted.bam//")
qualimap rnaseq -bam $file -gtf /data/$Annotation -outformat HTML -outdir /data/$QCAlignments/$DirName  --java-mem-size=5G
done')
$docker wait $DockerId


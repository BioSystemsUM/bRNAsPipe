#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline
# and creates a genome Index for use of STAR reads alignment tool
# it is called from script STAR.sh

### Variables
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='quay.io/biocontainers/star:2.7.3a--0'
BaseDir=$1
Genome=$2
Annotation=$3
Index=$4
ReadLength=$5
Threads=$6
((Param=ReadLength-1))
chmod 777 $BaseDir/$Index # give write permissions to directory that will have genome indexes
DockerId=$($docker run -d --rm --privileged -v $BaseDir:/data -u root -e Genome=$Genome -e Annotation=$Annotation -e Index=$Index -e Param=$Param -e Threads=$Threads $Image \
        STAR --runThreadN $Threads \
        --runMode genomeGenerate \
        --genomeFastaFiles /data/$Genome \
        --sjdbGTFfile /data/$Annotation \
        --sjdbOverhang $Param \
	--genomeDir /data/$Index)
$docker wait $DockerId

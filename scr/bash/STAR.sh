#!/usr/bin/env bash
# This script is part of RNA-seq analysis pipeline and it:
# - creates an STAR genome index if not already available
# - aligns reads

### Variables
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='quay.io/biocontainers/star:2.7.3a--0'
ScriptsFolder=$1
BaseDir=$2
Genome=$3
Annotation=$4
Index=$5
ReadLength=$6
ReadsDir=$7
AlignmentsDir=$8
SamplesInfoFile=$9
Study=${10}
Threads=${11}
SequencingType=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v study=$Study '$1==study {print $12}' | sort | uniq | head -1) # 'Unpaired' means single-end sequencing
                                                                                                                              # '1' means paired-end sequencing 

### Create a STAR genome index if not already available
if [[ ! -e $BaseDir/$Index ]] # if GenomeIndex directory doesn't exist we have to create it and create index
then
	mkdir $BaseDir/$Index
	$ScriptsFolder/STARGenomeIndex.sh $BaseDir $Genome $Annotation $Index $ReadLength $Threads
elif [[ -z `ls $BaseDir/$Index` ]] # if GenomeIndex directory is empty create index
then
        $ScriptsFolder/STARGenomeIndex.sh $BaseDir $Genome $Annotation $Index $ReadLength $Threads
fi

### Align reads to genome
# if study uses single-end sequencing:

if [[ $SequencingType == 'Unpaired' ]]
then
  Reads=$(cd $BaseDir/$ReadsDir ; ls | grep -E '.*\.fastq.gz$')
  for file in $Reads
  do Alignment=$(echo $file | sed 's/Trimmed.fastq.gz//')
     DockerId=$($docker run -u root --rm -d --privileged -v $BaseDir:/data $Image \
     STAR --runThreadN $Threads \
     --genomeDir /data/$Index \
     --readFilesIn /data/$ReadsDir/$file \
     --outFileNamePrefix /data/$AlignmentsDir/$Alignment \
     --readFilesCommand gunzip -c)
  $docker wait $DockerId
  done
else
  # if study uses paired-end sequencing:
  Samples=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v study=$Study '$1==study {print $11}' | tr -d '\r' | sort | uniq)
  for Sample in $Samples
  do echo $Sample
     fileF=$Sample"TrimmedFoward.fastq.gz"
     fileR=$Sample"TrimmedReverse.fastq.gz"
     echo $fileF
     echo $fileR
     DockerId=$($docker run -u root --rm -d --privileged -v $BaseDir:/data $Image \
     STAR --runThreadN $Threads \
     --genomeDir /data/$Index \
     --readFilesIn /data/$ReadsDir/$fileF /data/$ReadsDir/$fileR \
     --outFileNamePrefix /data/$AlignmentsDir/$Sample \
     --readFilesCommand gunzip -c)
  $docker wait $DockerId
  done
fi

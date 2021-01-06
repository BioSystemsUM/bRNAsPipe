#!/usr/bin/env bash
# This script downloads fastq raw files for all studies at same time
# and runs fastQC on raw data for those studies
# To Run this script: DownloadFiles.sh <study we want>
source ./Dirs.sh
ScriptsFolder=$BaseDir/'scr/bash'
SamplesInfoFile='data/Studies_RNAseq.txt'
Study=$1
StudyDir='data'/$Study
RawData='data'/$Study/'rawData'
FastqcRaw='data'/$Study/'fastQCRaw'
TrimOutDir='data'/$Study/'trimmedData'
FastqcTrim='data'/$Study/'fastQCtrim'
ReadsDir='data'/$Study/'trimmedData'
AlignmentsDir='data'/$Study/'Alignments'
CountsDir='data'/$Study/'RawCounts'  
NormDir='data'/$Study/'NormData'
# Create directories for a study in case they don't already exist:
source $ScriptsFolder/createDir.sh
DirArray=($BaseDir/$StudyDir \
          $BaseDir/$RawData \
          $BaseDir/$FastqcRaw \
          $BaseDir/$TrimOutDir \
          $BaseDir/$FastqcTrim \
          $BaseDir/$ReadsDir \
          $BaseDir/$AlignmentsDir \
          $BaseDir/$CountsDir \
	  $BaseDir/$NormDir)
createDirIfMissing ${DirArray[@]}

# Download raw data - fastq files:
cd $BaseDir/$StudyDir/'rawData' # change to directoy where we want to download raw data
cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v study=$Study '$1==study {print $13}' | tr -d '\r' | \
while read link
do wget -bqc $link 1>> $BaseDir/$StudyDir/$StudyId/'rawData/PIDs'
done

#!/usr/bin/env bash
# Run this script after downloading fastq raw data of one study, to determine read quality
# To Run this script: GetQCfiles.sh <study we want>
source Edit
ScriptsFolder=$BaseDir/'scr/bash'
Study=$1
SamplesInfoFile=$BaseDir/'data/StudiesRNAseq.txt'
RawData='data'/$Study/'rawData'
FastqcRaw='data'/$Study/'fastQCRaw'

$ScriptsFolder/fastqc.sh $BaseDir $RawData $FastqcRaw

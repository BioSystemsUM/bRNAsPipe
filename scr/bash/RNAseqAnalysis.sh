#!/usr/bin/env bash
# This script does RNAseq Analysis

### Variables
source ./Dirs.sh
ScriptsFolder=$BaseDir/'scr/bash'
Genome='GenomeRef'/$GenomeRef
Annotation='Annotation'/$AnnoFName
AnnotationDir=$(echo $Annotation | cut -d'/' -f 1)
Index='GenomeRef/Index'
SamplesInfoFile='data/Studies_RNAseq.txt'
Threads='40'
Study=$1
StudyDir='data'/$Study
RawData='data'/$Study/'rawData'
FastqcRaw='data'/$Study/'fastQCRaw'
TrimOutDir='data'/$Study/'trimmedData'
FastqcTrim='data'/$Study/'fastQCtrim'
ReadsDir='data'/$Study/'trimmedData'
AlignmentsDir='data'/$Study/'Alignments'
QCAlignments='data'/$Study/'QCAlignments'
CountsDir='data'/$Study/'RawCounts'
NormDir='data'/$Study/'NormData'
File4ReadLenCalc=$BaseDir/$RawData/$(ls $BaseDir/$RawData | grep -E ".*\.gz" | head -1)
source $ScriptsFolder/readLength.sh
ReadLength=$(ReadLenFct $File4ReadLenCalc)

### Run Trimmomatic to:
# - remove adapter sequences/primer sequences  and other contaminants (when they exist)
# - select reads with at least 36 bp length and 24 of average read quality score
$ScriptsFolder/trimmomatic.sh $BaseDir $RawData $TrimOutDir $Study $SamplesInfoFile $Threads
cd $BaseDir/$RawData # removes Raw fastq files to save space
rm $(find -name '*fastq.gz')

### Run FASTQC again - to evaluate read quality after running trimmomatic:
$ScriptsFolder/fastqc.sh $BaseDir $TrimOutDir $FastqcTrim

### Align reads with STAR:
$ScriptsFolder/STAR.sh $ScriptsFolder $BaseDir $Genome $Annotation $Index $ReadLength $ReadsDir $AlignmentsDir $SamplesInfoFile $Study $Threads
cd $BaseDir/$TrimOutDir # removes trimmed fastq files to save space
rm $(find -name '*fastq.gz')

### Run samtools to:
# - convert SAM files to BAM files (and then erase SAM files - it saves space)
# - sort BAM files by coordinates - sometimes needed for some downstream applications
# - create an index for BAM file - also sometimes needed for some downstream applications
$ScriptsFolder/samtools.sh $BaseDir $AlignmentsDir $Threads

### Evaluate alignments:
$ScriptsFolder/QCAlignments.sh $BaseDir $AlignmentsDir $QCAlignments $Annotation 

### Run htseq to count number of alignemnts to each gene:
$ScriptsFolder/htseq.sh $BaseDir $AlignmentsDir $Annotation $CountsDir

### Normalise counts of alinments to each gene:
$ScriptsFolder/NormCounts.sh $BaseDir $Annotation $AnnotationDir $SamplesInfoFile $Study $CountsDir $NormDir

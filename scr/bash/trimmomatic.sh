#!/usr/bin/env bash
# This script runs trimmomatic - RNAseqpipeline

### VARIABLES:
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='quay.io/biocontainers/trimmomatic:0.36--5'
BaseDir=$1
InDir=$2
OutDir=$3
Study=$4
SamplesInfoFile=$5
Threads=$6
SequencingType=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v study=$Study '$1==study {print $12}' | sort | uniq | head -1) # 'Unpaired' means single-end sequencing
                                                                                                                              # '1' means paired-end sequencing

### Trimmomatic Params - default:
Seq2RemoveParams='2:30:10'
LeadingParams='0'
TrailingParams='0'
SlidingWindowParams='0:0'
MinlenParams='36'
AvgQualParams='20'

# if study uses single-end sequencing:
if [[ $SequencingType == 'Unpaired' ]]
then
  for file in $(ls $BaseDir/$InDir | grep -E ".*\.gz")
  do
    OutFile=$(echo $file | sed 's/.fastq.gz/Trimmed.fastq.gz/')
    Seq2RemoveFile=$(echo $file | sed 's/.fastq.gz/Seq2Remove.fa/')
    TrimmomaticParamsFile=$(echo $file | sed 's/.fastq.gz/TrimParams/')
    OutSummary=$(echo $file | sed 's/.fastq.gz/Summary/')
    # if user defined sample trimmomatic params are available use them otherwise continue with default:
    [[ -f $BaseDir/$OutDir/$TrimmomaticParamsFile ]] && source $BaseDir/$OutDir/$TrimmomaticParamsFile
    # saves trimmomatic params (env variables) to a file that will need to use inside container (together with --env-file param):
    echo InDir=$InDir > Variables.txt
    echo OutDir=$OutDir >> Variables.txt
    echo Threads=$Threads >> Variables.txt
    echo Seq2RemoveParams=$Seq2RemoveParams >> Variables.txt
    echo LeadingParams=$LeadingParams >> Variables.txt
    echo TrailingParams=$TrailingParams >> Variables.txt
    echo SlidingWindowParams=$SlidingWindowParams >> Variables.txt
    echo MinlenParams=$MinlenParams >> Variables.txt
    echo AvgQualParams=$AvgQualParams >> Variables.txt
    echo file=$file >> Variables.txt
    echo OutFile=$OutFile >> Variables.txt
    echo Seq2RemoveFile=$Seq2RemoveFile >> Variables.txt
    echo OutSummary=$OutSummary >> Variables.txt
    # if file with contaminants/overrepresented sequences for that sample is available remove those sequences, otherwise don't remove sequences:
    [[ -f $BaseDir/$OutDir/$Seq2RemoveFile ]] && \
    DockerId=$($docker run --rm -d --privileged -u root --env-file Variables.txt -v $BaseDir:/data $Image bash -c 'trimmomatic SE -threads $Threads \
    /data/$InDir/$file /data/$OutDir/$OutFile \
    ILLUMINACLIP:/data/$OutDir/$Seq2RemoveFile:$Seq2RemoveParams \
    LEADING:$LeadingParams TRAILING:$TrailingParams SLIDINGWINDOW:$SlidingWindowParams MINLEN:$MinlenParams AVGQUAL:$AvgQualParams &> /data/$OutDir/$OutSummary') || \
    DockerId=$($docker run --rm -d --privileged -u root --env-file Variables.txt -v $BaseDir:/data $Image bash -c 'trimmomatic SE -threads $Threads \
    /data/$InDir/$file /data/$OutDir/$OutFile \
    LEADING:$LeadingParams TRAILING:$TrailingParams SLIDINGWINDOW:$SlidingWindowParams MINLEN:$MinlenParams AVGQUAL:$AvgQualParams &> /data/$OutDir/$OutSummary')
    $docker wait $DockerId
  done
else
# if study uses paired-end sequencing:
  Samples=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v study=$Study '$1==study {print $11}' | tr -d '\r' | sort | uniq)
  for Sample in $Samples
  do
    fileF=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v sample=$Sample '$11==sample {print $12,$13}' | awk '$1=="1" {print $2}' | awk -F'/' '{print $NF}')
    fileR=$(cat $BaseDir/$SamplesInfoFile | awk -F '\t' -v sample=$Sample '$11==sample {print $12,$13}' | awk '$1=="2" {print $2}' | awk -F'/' '{print $NF}')
    OutFileFpaired=$(echo $Sample"TrimmedFoward.fastq.gz")
    OutFileFunpaired=$(echo $Sample"TrimmedFowardUnpaired.fastq.gz")
    OutFileRpaired=$(echo $Sample"TrimmedReverse.fastq.gz")
    OutFileRunpaired=$(echo $Sample"TrimmedReverseUnpaired.fastq.gz")
    Seq2RemoveFile=$(echo $Sample"Seq2Remove.fa")
    TrimmomaticParamsFile=$(echo $Sample"TrimParams")
    OutSummary=$(echo $Sample"Summary")
    # if user defined sample trimmomatic params are available use them othewise continue with default:
    [[ -f $BaseDir/$OutDir/$TrimmomaticParamsFile ]] && source $BaseDir/$OutDir/$TrimmomaticParamsFile
    # saves trimmomatic params (env variables) to a file that will need to use inside container (together with --env-file param):
    echo InDir=$InDir > Variables.txt
    echo OutDir=$OutDir >> Variables.txt
    echo Threads=$Threads >> Variables.txt
    echo Seq2RemoveParams=$Seq2RemoveParams >> Variables.txt
    echo LeadingParams=$LeadingParams >> Variables.txt
    echo TrailingParams=$TrailingParams >> Variables.txt
    echo SlidingWindowParams=$SlidingWindowParams >> Variables.txt
    echo MinlenParams=$MinlenParams >> Variables.txt
    echo AvgQualParams=$AvgQualParams >> Variables.txt
    echo fileF=$fileF >> Variables.txt
    echo fileR=$fileR >> Variables.txt
    echo OutFileFpaired=$OutFileFpaired >> Variables.txt
    echo OutFileFunpaired=$OutFileFunpaired >> Variables.txt
    echo OutFileRpaired=$OutFileRpaired >> Variables.txt
    echo OutFileRunpaired=$OutFileRunpaired >> Variables.txt
    echo Seq2RemoveFile=$Seq2RemoveFile >> Variables.txt
    echo OutSummary=$OutSummary >> Variables.txt
    #if file with contaminants/overrepresented sequences for that sample is available remove those sequences, otherwise don't remove sequences:
    [[ -f $BaseDir/$OutDir/$Seq2RemoveFile ]] && \
    DockerId=$($docker run --rm -d --privileged -u root --env-file Variables.txt -v $BaseDir:/data $Image bash -c 'trimmomatic PE -threads $Threads \
    /data/$InDir/$fileF /data/$InDir/$fileR /data/$OutDir/$OutFileFpaired /data/$OutDir/$OutFileFunpaired /data/$OutDir/$OutFileRpaired /data/$OutDir/$OutFileRunpaired \
    ILLUMINACLIP:/data/$OutDir/$Seq2RemoveFile:$Seq2RemoveParams:2:keepBothReads \
    LEADING:$LeadingParams TRAILING:$TrailingParams SLIDINGWINDOW:$SlidingWindowParams MINLEN:$MinlenParams AVGQUAL:$AvgQualParams &> /data/$OutDir/$OutSummary') || \
    DockerId=$($docker run --rm -d --privileged -u root --env-file Variables.txt -v $BaseDir:/data $Image bash -c 'trimmomatic PE -threads $Threads \
    /data/$InDir/$fileF /data/$InDir/$fileR /data/$OutDir/$OutFileFpaired /data/$OutDir/$OutFileFunpaired /data/$OutDir/$OutFileRpaired /data/$OutDir/$OutFileRunpaired \
    LEADING:$LeadingParams TRAILING:$TrailingParams SLIDINGWINDOW:$SlidingWindowParams MINLEN:$MinlenParams AVGQUAL:$AvgQualParams &> /data/$OutDir/$OutSummary')
    $docker wait $DockerId
  done
fi
rm Variables.txt

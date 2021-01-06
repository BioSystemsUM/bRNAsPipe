#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline
# It normalizes raw counts of alignments to each gene

### Variables:
docker='/usr/bin/podman' # turing works with podman instead of docker
ImageGTF='tlbbdocker/gtftools:v1'
ImageNorm='tlbbdocker/normcounts:v1'
BaseDir=$1
Annotation=$2
AnnotationDir=$3
SamplesInfoFile=$4
Study=$5
CountsDir=$6
NormDir=$7
AnnotationTmp=$(echo $Annotation"Tmp")

### Get file with gene length from gtf file:
# Copy gtf, replace 1st field by value '1' and use tabs to separate fields (explanation below):
grep -v "^#" $BaseDir/$Annotation | awk '{OFS="\t"} $1=1' > $BaseDir/$AnnotationTmp
# Calculate gene length:
DockerIdGTF=$($docker run --privileged -d -u root --rm -v $BaseDir:/data -e AnnotationTmp=$AnnotationTmp -e AnnotationDir=$AnnotationDir $ImageGTF \
bash -c "./gtftools.py -l /data/$AnnotationDir/GeneLengthTmp /data/$AnnotationTmp ; \
cat /data/$AnnotationDir/GeneLengthTmp | cut -f1,5 | sed 's/merged/length/' > /data/$AnnotationDir/GeneLength.bed ; \
rm /data/$AnnotationDir/GeneLengthTmp ; \
rm /data/$AnnotationTmp")
$docker wait $DockerIdGTF
# Note: - gtftools has a hitch. If chromossome names in gtf file is not on a specific format, does not calculate length of genes on non-autossomic chromossomes.
#         So, we have to exclude commented lines and change the chromossome field (1st field) from chrN (in case of NCBI gtf) to '1'
#       - length is also not determined for pseudogenes in case of gtf NCBI files 
#       - see this link: http://genomespot.blogspot.com/2019/01/using-gtf-tools-to-get-gene-lengths.html

### Join datasets' raw counts into a table (with samples of same condition side by side) and normalize:
# Save env variables/params to a file so that we can use them inside a container (together with --env-file):
echo Study=$Study > Variables
echo SamplesInfoFile=$SamplesInfoFile >> Variables
echo CountsDir=$CountsDir >> Variables
echo AnnotationDir=$AnnotationDir >> Variables
echo NormDir=$NormDir >> Variables

DockerIdNorm=$($docker run --privileged -d --rm --env-file Variables -v $BaseDir:/data $ImageNorm \
bash -c 'Rscript /data/scr/R/NormCountsRNAseq.R \
/data/scr/R/NormCountsRNAseqFunctions.R \
$Study \
/data \
$SamplesInfoFile \
$CountsDir \
$AnnotationDir \
$NormDir')
$docker wait $DockerIdNorm
rm Variables

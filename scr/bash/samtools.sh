#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline and:
# - converts SAM files to BAM files (and then erases SAM files - it saves space)
# - sorts BAM files by coordinates - sometimes needed for some downstream applications
# - creates an index for BAM file - also sometimes needed for some downstream applications

### Variables:
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='biocontainers/samtools:v1.3_cv2'
BaseDir=$1
AlignmentsDir=$2
Threads=$3
DockerId=$($docker run --privileged --rm -d -u root -e Threads=$Threads -v $BaseDir/$AlignmentsDir:/data $Image \
bash -c 'cd /data ; \
for file in $(ls | grep -E '.*Aligned.out.sam$') ; \
do fileOutSort=$(echo $file | sed 's/Aligned.out.sam/.Sorted.bam/') ; \
  samtools view -bS /data/$file | samtools sort - -@$Threads -m 10G -o /data/$fileOutSort ; \
  samtools index /data/$fileOutSort ; \
done ; \
rm *.sam')
$docker wait $DockerId

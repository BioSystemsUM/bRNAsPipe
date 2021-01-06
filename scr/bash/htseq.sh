#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline and runs Htseq to count number of reads aligning to each gene:

### Variables
docker='/usr/bin/podman' # turing works with podman instead of docker
Image='biocontainers/htseq:v0.11.2-1-deb-py3_cv1'
BaseDir=$1
AlignmentsDir=$2
Annotation=$3
CountsDir=$4
DockerId=$($docker run --rm --privileged -d -u root -e AlignmentsDir=$AlignmentsDir -e Annotation=$Annotation -e CountsDir=$CountsDir -v $BaseDir:/data $Image \
bash -c 'cd /data/$AlignmentsDir ; \
for file in $(ls | grep -E ".*Sorted.bam$") ; \
do Count=$(echo $file | sed "s/Sorted.bam/RawCounts/") ; \
   Name=$(echo $Count | cut -d"." -f1) ; \
   htseq-count -f bam -m=intersection-strict --stranded=no -t gene -i gene_id /data/$AlignmentsDir/$file /data/$Annotation > /data/$CountsDir/$Count"Tmp" ; \
   cat /data/$CountsDir/$Count"Tmp" | grep -E "__no_feature|__ambiguous|__too_low_aQual|__not_aligned|__alignment_not_unique" > /data/$CountsDir/$Name"Summary" ; \
   cat /data/$CountsDir/$Count"Tmp" | grep -Ev "__no_feature|__ambiguous|__too_low_aQual|__not_aligned|__alignment_not_unique" > /data/$CountsDir/$Count ; \
done; \
cd /data/$CountsDir ; rm *Tmp')
$docker wait $DockerId

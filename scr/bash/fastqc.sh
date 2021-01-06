#!/usr/bin/env bash
# This script is part of RNAseq analysis pipeline and applies FASTQC 
docker='/usr/bin/podman' # turing works with podman instead of docker
Image="biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"
BaseDir=$1
In=$2
Out=$3
DockerId=$($docker run -d --rm -u root --privileged -v $BaseDir:/data -e In=$In -e Out=$Out $Image bash -c 'cd /data/$In ; for file in $(ls | grep -E ".*\.gz") ; do fastqc $file &> /dev/null ; done ;  mv $(ls | grep -E ".html|.zip") /data/$Out')
$docker wait $DockerId

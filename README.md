## Pipeline for RNAseq and Microarray analysis
### RNAseq
Pipeline for RNAseq was developed in Bash and it uses docker containers. Requirements to run are: *Linux* system and *Podman*. It is recommended to use ensembl annotation and genome refference.
 - Example annotation file: [Homo_sapiens.GRCh38.99.gtf.gz](ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz)
 - Example genome reference fasta: [Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
1. To download files above, move to Annotation and GenomeRef folders respectively and do: `wget <link>`
2. To unzip files above: `gunzip <FileName>`
3. Define base directory and genome and annotation files path in **scr/bash/Dirs.sh**
4. Table **Studies_RNAseq.txt** is a tab-delimited file, under **data** directory, that should be filled with studies info.
 - *Study* is the study identifier
 - *SampleId* is sample identifier
 - *Reads* has '1' or '2' to destinguish between foward and reverse reads, single-end studies have 'Unpaired'
 - *Link* is the link to fastq.gz file.
 - All other columns should be filled with 'NA' when there are no values.
5. Download files of a study with:  `DownloadFiles.sh <Study>`
6. Confirm if files finished to download: `ps -e | grep <jobId>` To get Job ids of donwloads `cd /data/<Study>/rawData` and do `cat PIDs` 
7. Get fastQC files for a study with: `GetQCfiles.sh <Study>`
8. After this, manually check fastqc results and decide which contaminants/overrepresented sequences should be removed in each sample and add them to file **<SampleName>Seq2RemoveFile** in folder **trimmedData** so that Trimmomatic will remove those sequences. If no file is provided, trimmomatic runs without excluding those sequences. Example of **<SampleName>Seq2RemoveFile** content:
>seqname
ACTTTTTTTTTTTTTTTTTTT 
9. To define specific trimmomatic parameters for a sample, include a file named **<sampleName>TrimParams** in directory **trimmedData** where you can change trimmomatic parameters for each sample, if you see for example that reads need to be trimmed in that study. Otherwise, default parameters are run.
10. To run the rest of the analysis: `RNAseqAnalysis.sh <Study>`
### Microarray
To run in Windows OS with R.
File with studies info is: *Studies_Microarrays.xlsx*
Run script *scr/R/MicroarrayNormalize.R*
Paths are hardcoded.


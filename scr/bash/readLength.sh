#!/usr/bin/env bash
# This script has a function that determines gene length of a fastq.gz file
# Sometimes read length varies even within same sample
# so, we make a distribution of read length in a file a choose the length that is shared by more reads
ReadLenFct ()
{
    FastqGzFile=$1
    distLen=$(gunzip -c $FastqGzFile | awk '{if(NR%4==2) print length($1)}' | sort | uniq -c ) # NR%4 numbers lines with repetitions of numbers 0,1,2,3
                                                                                               # if(NR%4==2) selects every second line (which is numbered with 2) which corresponds to a read
	    										       # print length($1) prints the lenth of 1st field - read length
											       # sort sorts read lengths from lower value to higher values
											       # uniq -c counts number of times each read length exists in a file
    FrqReadLen=$(echo $distLen | sort -nr -s -k 1,1 | head -1 | cut -d' ' -f2) # -nr sorts by number in reverse order and -s says sort just by the specified field (which is defined by -k 1,1)
                                                                               # this gives the all line sorted just by first field (frequency of read) fom high to low value
               				                                       # head -1 selects the line with most frequent read length
				                                               # cut -d ' ' -f2 selects the most frequent read length
    echo $FrqReadLen
}

#!/bin/sh

##### This script is used for the detection of uniquely-mapped paired-end viral reads
#####
#
# Author:       Xun Chen
# Email:        xunchen85@gmail.com
# Date:         09/10/2020

##### Required dependencies in the enviroment
# bwa-0.7.17
# samtools-1.10 
# Repeatmasker latest version
# trf latest version
# dust latest version

##### Input file
# Clean paired-end unmapped reads with naming of sampleID_unmapped_1.fq and sampleID_unmapped_2.fq

##### Variables
sampleID=$1
threads=$2

##### Required databases
virus_genome="virus_db_920k.fa"		### (should be indexed by using bwa)
virus_size="virus_db_920k.virus_size"
virus_list="virus_db_920k.virus_list"
vector_genome="vector_db.fa"		### (should be indexed by using bwa)
human_genome="human_genome.fa"		### (should be indexed by using bwa)

######### Section 1 Detect viral reads for each sample
##### Step 1 Align the paired-end unmapped reads against the viral genome
bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $virus_genome ${sampleID}_unmapped_1.fq ${sampleID}_unmapped_2.fq >${sampleID}_virus.sam
perl Reformat_virusSam1.pl ${sampleID}_virus.sam 30 >${input_sampleID}.read1
perl Reformat_virusSam2.pl ${input_sampleID}.read1 $virus_list >${input_sampleID}.reads2

##### Step 2 Re-align the paired-end unmapped reads against the human genome
bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $human_genome ${sampleID}_unmapped_1.fq ${sampleID}_unmapped_2.fq >${sampleID}_human.sam
perl Remove_mappedReads.pl ${sampleID}_human.sam | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk '{if ($5-$4+1>=50 || $8-$7+1>=50)print$0}' |awk '{print$2}' | sort | uniq > ${sampleID}.human

##### Step 3 Align the paired-end unmapped reads against the vector database (AS >= 50)
##### Check vector sequence
bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $vector_genome ${sampleID}_unmapped_1.fq ${sampleID}_unmapped_2.fq >${sampleID}_vector.sam
perl Remove_mappedReads.pl ${sampleID}_vector.sam | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk '{if ($5-$4+1>=50 || $8-$7+1>=50)print$0}' |awk '{print$2}' | sort | uniq > ${sampleID}.vector

##### Step 4 Screen for repeats
perl Convert.pl ${sampleID}
trf ${input2}.fa 2 5 7 80 10 50 2000 -h -f -d -m
dust ${input2}.fa.2.5.7.80.10.50.2000.mask >${input2}.fa2
perl Remove_reads_trf_dust.pl ${input2}.fa2 >${input2}.fa3
RepeatMasker -e rmblast -pa $threads -s -species human ${sampleID}.fa3
perl Reformat_repeatmaskerOutput.pl ${sampleID} >${sampleID}.repeat
perl Reformat_repeats.pl ${sampleID} >${sampleID}.repeat2

##### Step 5 Filter and Organize the output file
perl Filtering_false_reads.pl ${sampleID} vhr >${input_sampleID}.virus
perl Filtering_virus_v.pl ${sampleID}.virus $virus_size >${sampleID}.virus_v
perl Filtering_virus_g.pl ${sampleID}.virus $virus_size >${sampleID}.virus_g
perl Summary_virus1.pl ${sampleID}.virus_v ${sampleID}.virus_g $virus_size $virus_list >${sampleID}.virus_v2
perl Summary_virus2.pl ${sampleID}.virus_v2 0 | uniq >${sampleID}.virusAbundance
rm ${sampleID}.virus ${sampleID}.virus_v ${sampleID}.virus_g ${sampleID}.virus_v2

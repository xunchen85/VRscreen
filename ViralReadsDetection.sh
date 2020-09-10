#!/bin/sh

##### This script is used for the detection of uniquely-mapped viral reads
#####
#
# Author:       Xun Chen
# Email:        xunchen85@gmail.com
# Date:         09/10/2020

##### Required dependencies in the enviroment
bwa
samtools 
repeatmasker

##### Required databases
virus.fa
virus.fa.virus_list
virus_size
human_genome.fa
vector.fa

##### Input file
Clean paired-end unmapped reads

##### Variables
threads=$1
sampleID=$2

##### Section 1
bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $virus_genome ${sampleID}_unmapped_1.fq ${sampleID}_unmapped_2.fq >${sampleID}_vsu.sam

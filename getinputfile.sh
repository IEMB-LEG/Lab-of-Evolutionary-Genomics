#! /bin/perl -w
open INFO, "list" or die "cannot open file: $!"; # list file contains line fq names, without the _fq
while ($i=<INFO>) {
chomp ($i);
open OUT, ">$i.mpileup.bwa" or die "cannot open file: $!";
print OUT '#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=4,vmem=15g,walltime=48:00:00
#PBS -M 1534732467@qq.com
#PBS -m abe
#PBS -N '.$i."\n".'#PBS -j oe'."\n".'cd /N/dc2/projects/MicroEukMA/students/eric/Shewanella/'."\n".'/N/dc2/projects/MicroEukMA/softwares/samtools-1.3.1/samtools mpileup -s -f ../she.fasta '.$i.'.sort.dedup.bam > '.$i.'.dedup.mpileupBwa'."\n".'perl /N/dc2/projects/MicroEukMA/softwares/my_working_scripts/printBasesBwa.pl '.$i.'.dedup.mpileupBwa > '.$i.'.printBasesBwa';

#! /bin/bash
#PBS -q batch
#PBS -k o
#PBS -l nodes=1:ppn=1,vmem=50g,walltime=10:00:00
#PBS -M wukun@stu.ouc.edu.cn
#PBS -m abe
#PBS -N MMR_coverage
#PBS -j oe
cd /N/dc2/projects/MicroEukMA/students/eric_1/eric/Shewanella/
R CMD BATCH sel-linux.R

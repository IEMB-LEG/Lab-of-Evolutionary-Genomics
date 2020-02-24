#! /bin/perl -w
# do not do too many samples at one time, do 40 samples batches
$file=$ARGV[0];
$name=$file;
open INFO, "$file" or die "cannot open file: $!"; # list file contains line fq names, without the _fq
# to avoid java outof memory in MarkDuplicate, do: -xmx2g or-XX:-UseGCOverheadLimit
@array=();
while ($i=<INFO>) {
chomp ($i);
push @array, "--variant $i ";
}
print "#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -M lc9869\@stu.ouc.edu.cn
#PBS -m abe
#PBS -N ".$name."
#PBS -j oe
cd /N/dc2/projects/MicroEukMA/students/panjiao/"."\n";
print "module load java/jre/1.8.0_73"."\n";
print "/N/dc2/projects/MicroEukMA/softwares/samtools-1.3.1/samtools faidx photo.fasta"."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R= photo.fasta O= photo.dict"."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -T CombineGVCFs -R photo.fasta ";
foreach (@array) {
print "$_";
}
print "-o ".$name.".merge.g.vcf"."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R photo.fasta --variant ".$name.".merge.g.vcf -o ".$name.".raw.vcf "."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -T SelectVariants -R photo.fasta -V ".$name.".raw.vcf -selectType SNP -o ".$name.".snp.vcf"."\n";
print 'java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -R photo.fasta -T VariantFiltration -V '.$name.'.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o '.$name.'.filtered_snp.vcf'."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -T SelectVariants -R photo.fasta -V ".$name.".raw.vcf -selectType INDEL -o ".$name.".indel.vcf"."\n";
print 'java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -T VariantFiltration -R photo.fasta -V '.$name.'.indel.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o '.$name.'.filtered_indel.vcf'."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -R photo.fasta -T VariantsToTable -V ".$name.".filtered_snp.vcf -F CHROM -F POS -F REF -F ALT -F HOM-REF -F HOM-VAR -o ".$name."_snp.table"."\n";
print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -R photo.fasta -T VariantsToTable -V ".$name.".filtered_indel.vcf -F CHROM -F POS -F REF -F ALT -F HOM-REF -F HOM-VAR -o ".$name."_indel.table"."\n";
#print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -R photo.fasta -T VariantsToTable -V ".$name.".filtered_snp.vcf -F CHROM -F POS -F REF -F ALT -F HOM-REF -F HOM-VAR -o ".$name."_snp.table"."\n";
#print "java -Xmx2g -jar /N/dc2/projects/MicroEukMA/softwares/GATK_3.6/GenomeAnalysisTK.jar -R photo.fasta -T VariantsToTable -V ".$name.".filtered_indel.vcf -F CHROM -F POS -F REF -F ALT -F HOM-REF -F HOM-VAR -o ".$name."_indel.table";

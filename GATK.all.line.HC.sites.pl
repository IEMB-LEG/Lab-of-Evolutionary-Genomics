#! /bin/env perl
# file is all.uni.vcf
use strict;
use warnings;
my $file=$ARGV[0];
open(FILE,"$file") || die "cannot open file";
my @line=();
my %hash_site=();
my %hash_A=();
my %hash_C=();
my %hash_G=();
my %hash_T=();
my %hash_coverage=();

while (<FILE>) {
	chomp $_;
	my @input=split(/\s+/,$_);
	if ($input[0] =~ "##") { next;}
	if ($input[0] =~ "CHROM") { @line=@input[9..$#input];
					next;}
#	if ($input[7] eq '.') {next;}
#	my @info=split(/;/,$input[7]);
#	my @MQ_info=split(/=/,$info[2]);
#	my $mq=$MQ_info[1];
	my $counter=-1;
#	if($input[5] eq ".") {next;} ;
#	if (($input[5] <100) or ($mq < 59) or ()) {next;}
	if ($input[6] ne "PASS") {next;} 	
	my @reference=split(//,$input[3]);
	foreach my $sample (@input[9..$#input]) {
			$counter++;

#			if (($sample =~ ":") and ($sample !~ ",")) {
#				my @COV=split(/:/,$sample);
#				my $coverage=$COV[1];
#				$hash_site{$line[$counter]}+=1;
#				if ($coverage eq ".") { $coverage=0;}
#				$hash_coverage{$line[$counter]}+=$coverage;
#				}
#			if (($sample =~ ":") and ($sample =~ ",")) {
                                my @elements=split(/:/,$sample);
                                my $coverage=$elements[2];
				
				if ($sample !~ /\.\/\./) {
                                $hash_site{$line[$counter]}+=1;
                                $hash_coverage{$line[$counter]}+=$coverage;
                               
			if ($reference[0] =~ "A") {$hash_A{$line[$counter]}+=1;}
			if ($reference[0] =~ "C") {$hash_C{$line[$counter]}+=1;}
			if ($reference[0] =~ "G") {$hash_G{$line[$counter]}+=1;}
			if ($reference[0] =~ "T") {$hash_T{$line[$counter]}+=1;}
			}
	}
}	

while( my ($k, $v) = each %hash_site ) {
        print "Line: $k Coverage: ".$hash_coverage{$k}/$v." A_sites: ".$hash_A{$k}." C_sites: ".$hash_C{$k}." G_sites: ".$hash_G{$k}." T_sites: ".$hash_T{$k}."\n";
}

#! /bin/env perl
# curate out heavily affected by mismappin sites first
use strict;
use warnings;
my $file=$ARGV[0];
open(FILE,"$file") || die "cannot open file";
my @line=();
while (<FILE>) {
	chomp $_;
	my @input=split(/\s+/,$_);
	if ($input[0] =~ "#CHROM") { @line=@input[3..$#input];
					next;}
	my $counter=-1;
	foreach my $sample (@input[3..$#input]) {
			$counter++;	
			if (($sample ne $input[2]) and ($sample ne "-")) {print $line[$counter]." ".$input[0]." ".$input[1]." ".$input[2]." ".$sample."\n";}
				}
}


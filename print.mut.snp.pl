#! /bin/env perl
use strict;
use warnings;
my $file=$ARGV[0];
open(FILE,"$file") || die "cannot open file";
while (<FILE>) {
	chomp;
	my @input=split(/\s+/,$_);
	my $mut=$input[0]+4;
	print $input[0]." ".$input[1]." ".$input[2]." ".$input[$#input]." ".$input[$mut]."\n";
}

